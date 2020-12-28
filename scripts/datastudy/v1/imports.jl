ENV["GKSwstype"] = "nul"  # For StatsPlots

import Pkg; Pkg.activate("../../../")
include("../../info.jl")

using CytofDiffDensity; const cdd = CytofDiffDensity
using DataFrames, CSV
using Util
import Random
using DrWatson
using MCMC
using BSON
using LaTeXStrings
using StatsPlots
using StatsFuns

# required for `BSON.load`
using Distributions, MCMC
import MCMC: Bijectors

# Plot settings
plotsize = (450, 450)
Plots.scalefontsizes()
Plots.scalefontsizes(1.5)

simname = "v1"
make_resultsdir(sim) = mkpath("$(Info.resultsdir_datastudy)/$(simname)/$(savename(sim))")
data_path = "../../../data/TGFBR2/cytof-data/donor1.csv"

read_data(marker::Symbol) = read_data(data_path, marker)
function read_data(path_to_data::String, marker::Symbol;
                   treatmentsym::Symbol=:treatment, log_response::Bool=true)
  data = DataFrame(CSV.File(path_to_data))
  yC_all = data[ismissing.(data[!, treatmentsym]), marker]
  yT_all = data[.!ismissing.(data[!, treatmentsym]), marker]

  if log_response
    yC_all .= log.(yC_all)
    yT_all .= log.(yT_all)
  end

  yC = yC_all[isfinite.(yC_all)]
  yT = yT_all[isfinite.(yT_all)]

  NC = length(yC_all) 
  NT = length(yT_all) 
  QC = NC - length(yC)
  QT = NT - length(yT)

  return (yC=yC, yT=yT, NC=NC, NT=NT, QC=QC, QT=QT)
end

function plot_gtilde_data!(data; ygrid=nothing, colorC=:blue, colorT=:red)
  histogram!(data.yC, normalize=true, label=nothing, la=0, alpha=.3, color=colorC)
  histogram!(data.yT, normalize=true, label=nothing, la=0, alpha=.3, color=colorT)
end

function run(sim)
  sim = (; sim...)

  # Simulate data.
  data = read_data(sim.marker)

  # Make results dir if needed.
  resultsdir = make_resultsdir(sim)

  # Create Priors.
  y = [data.yC; data.yT]

  # Plot and save data
  imgdir = mkpath(joinpath(resultsdir, "img"))
  plot(size=plotsize)
  plot_gtilde_data!(data)
  savefig(joinpath(imgdir, "data.pdf"))
  closeall()

  Random.seed!(0)
  mu_prior = Normal(mean(y), std(y))  # independent priors for mu_k.
  if sim.skewtmix
    # SkewT mixture
    model = Gtilde(data.yC, data.yT, sim.K, true, mu=mu_prior, skew=true, tdist=true)
  else
    # Normal mixture
    model = Gtilde(data.yC, data.yT, sim.K, true, mu=mu_prior, skew=false, tdist=false)
  end
  exclude = [:vC, :vT, :zetaC, :zetaT, :lambdaC, :lambdaT]
                 
  init = MCMC.make_init_state(model)
  spl = make_sampler(model, init=init)

  nburn, nsamps, thin = (30000, 4000, 2)
  # nburn, nsamps, thin = (3, 4, 2)
  callback = make_callback(model, nburn=nburn, nsamps=nsamps, thin=thin)

  Util.redirect_stdout_to_file(joinpath(resultsdir, "log.txt")) do
    # Run MCMC.
    chain, metrics, laststate = mcmc(spl, nsamps, init=init, nburn=nburn,
                                     thin=thin, callback=callback,
                                     exclude=exclude, return_last_state=true);

    # Save chain and metrics.
    BSON.bson(joinpath(resultsdir, "results.bson"),
              Dict(:chain => chain, :metrics => metrics, :data => data,
                   :model => model, :laststate => laststate, :spl => spl))
  end
end

function _postprocess(sim, resultsdir)
  sim = (; sim...)

  # Make image dir if needed.
  imgdir = mkpath(joinpath(resultsdir, "img"))

  # Load results.
  # chain, metrics, data, model.
  r = (; BSON.load(joinpath(resultsdir, "results.bson"))...)

  # Print summary stats
  cdd.printsummary(r.chain, r.metrics); println()

  # Print data info
  foreach(s -> println(s, ": ", getindex(r.data, s)), [:QC, :QT, :NC, :NT])

  # Print model info
  cdd.print_model_info(r.model)

  # Plot loglike
  plot(r.metrics[:loglike], label=nothing)
  savefig(joinpath(imgdir, "loglike.pdf"))
  closeall()

  # Plot posterior for each multivariate parameter.
  for sym in [:mu, :sigma, :nu, :phi, :etaC, :etaT]
    p = hcat(getindex.(r.chain, sym)...)'

    plot(p, label=nothing)
    savefig(joinpath(imgdir, "$(sym)-trace.pdf"))
    closeall()

    boxplot(p, label=nothing)
    savefig(joinpath(imgdir, "$(sym).pdf"))
    closeall()
  end

  # Plot posterior for each univariate parameter.
  for sym in [:tau]
    p = vec(hcat(getindex.(r.chain, sym)...))

    plot(p, label=nothing)
    savefig(joinpath(imgdir, "$(sym)-trace.pdf"))
    closeall()

    histogram(p, label=nothing, normalize=true)
    savefig(joinpath(imgdir, "$(sym).pdf"))
    closeall()
  end

  # Plot posterior density for C and T.
  y = [r.data.yC; r.data.yT]
  ygrid = cdd.make_ygrid(y, 100)
  plot(size=plotsize)
  cdd.plot_post_density!(r.chain, ygrid)
  plot_gtilde_data!(r.data, colorC=:grey, colorT=:grey)
  savefig(joinpath(imgdir, "post-density.pdf"))
  closeall()

  # Plot gammas (zero part)
  begin
    pzero = Pzero(NC=r.data.NC, NT=r.data.NT,
                  QC=r.data.QC, QT=r.data.QT)
    pzero_inference = cdd.infer_Pzero1(pzero, 10000)
    gammaC_post = pzero_inference.distC
    gammaT_post = pzero_inference.distT
    gamma_grid = cdd.make_gamma_grid(gammaC_post, gammaT_post, a=.001)

    plot(size=plotsize)
    plot!(gamma_grid, pdf.(gammaC_post, gamma_grid), label=nothing, color=:blue)
    plot!(gamma_grid, pdf.(gammaT_post, gamma_grid), label=nothing, color=:red)
    cdd.plot_gamma_uq!(gammaC_post, color=:blue, Q=pzero.QC, N=pzero.NC)
    cdd.plot_gamma_uq!(gammaT_post, color=:red, Q=pzero.QT, N=pzero.NT)
    savefig(joinpath(imgdir, "gamma.pdf"))
    closeall()
  end

  # Plot F̃ᵢ
  begin 
    # Plot Fi-tilde
    plot(size=plotsize)
    cdd.plot_Fi_tilde_cdf!(r.model, r.chain, pzero)
    ylims!(0, 1)
    xlabel!(L"\tilde{y}")
    ylabel!("CDF")
    savefig(joinpath(imgdir, "Fi-tilde-postmean.pdf"))
    closeall()

    # Plot Fi
    plot(size=plotsize)
    cdd.plot_Fi_tilde_cdf!(r.model, r.chain, pzero, exponentiate=true)
    ylims!(0, 1)
    xlabel!("y")
    ylabel!("CDF")
    savefig(joinpath(imgdir, "Fi-postmean.pdf"))
    closeall()
  end
end


function postprocess(sim)
  # Make results dir if needed.
  resultsdir = make_resultsdir(sim)

  Util.redirect_stdout_to_file(joinpath(resultsdir, "info.txt")) do
    _postprocess(sim, resultsdir)
  end
end


function parse_dic(model_info::String)
  m = match(r"(?<=DIC:\s)\d+\.\d+", model_info).match
  return parse(Float64, m)
end


# TODO:
# - [ ] Implement `resume`.
