# TODO:
# - [ ] use independent priors for mu (no ordering)

ENV["GKSwstype"] = "nul"  # For StatsPlots

import Pkg; Pkg.activate("../../../")
include("../../info.jl")
include("scenarios.jl")

using CytofDiffDensity; const cdd = CytofDiffDensity
using Util
import Random
using DrWatson
using MCMC
using BSON
using LaTeXStrings
using Distributions  # required for `BSON.load`
using StatsPlots
using StatsFuns

# Plot settings
plotsize = (450, 450)
Plots.scalefontsizes()
Plots.scalefontsizes(1.5)

simname = "v3"
make_resultsdir(sim) = mkpath("$(Info.resultsdir_simstudy)/$(simname)/$(savename(sim))")

function plot_simdata!(simdata; ygrid=nothing)
  histogram!(simdata.yC, normalize=true, label=nothing, la=0, alpha=.3, color=:blue)
  histogram!(simdata.yT, normalize=true, label=nothing, la=0, alpha=.3, color=:red)
  ygrid === nothing && (ygrid = cdd.make_ygrid([simdata.yC; simdata.yT], 100))
  plot!(ygrid, pdf.(simdata.mmC, ygrid), label=nothing, color=:blue)
  plot!(ygrid, pdf.(simdata.mmT, ygrid), label=nothing, color=:red)
end

function run(sim)
  sim = (; sim...)

  # Simulate data.
  simdata = scenarios(sim.snum, Ni=100_000, seed=1)

  # Make results dir if needed.
  resultsdir = make_resultsdir(sim)

  # Create Priors.
  y = [simdata.yC; simdata.yT]

  # Plot and save data
  imgdir = mkpath(joinpath(resultsdir, "img"))
  plot(size=plotsize)
  plot_simdata!(simdata)
  savefig(joinpath(imgdir, "data.pdf"))
  closeall()

  Random.seed!(0)
  mu_prior = Normal(mean(y), std(y))  # independent priors for mu_k.
  if sim.skewtmix
    # SkewT mixture
    model = Gtilde(simdata.yC, simdata.yT, sim.K, true, mu=mu_prior, skew=true, tdist=true)
  else
    # Normal mixture
    model = Gtilde(simdata.yC, simdata.yT, sim.K, true, mu=mu_prior, skew=false, tdist=false)
  end
  exclude = [:vC, :vT, :zetaC, :zetaT, :lambdaC, :lambdaT]
                 
  init = MCMC.make_init_state(model)
  spl = make_sampler(model, init=init)

  nburn, nsamps, thin = (15000, 4000, 2)
  callback = make_callback(model, nburn=nburn, nsamps=nsamps, thin=thin)

  Util.redirect_stdout_to_file(joinpath(resultsdir, "log.txt")) do
    # Run MCMC.
    chain, metrics = mcmc(spl, nsamps, init=init, nburn=nburn, thin=thin,
                          callback=callback, exclude=exclude);

    # Save chain and metrics.
    BSON.bson(joinpath(resultsdir, "results.bson"),
              Dict(:chain => chain, :metrics => metrics, :simdata => simdata,
                   :model => model))
  end
end

function _postprocess(sim, resultsdir)
  sim = (; sim...)

  # Make image dir if needed.
  imgdir = mkpath(joinpath(resultsdir, "img"))

  # Load results.
  # chain, metrics, simdata, model.
  r = (; BSON.load(joinpath(resultsdir, "results.bson"))...)

  # Print summary stats
  cdd.printsummary(r.chain, r.metrics); println()

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
  y = [r.simdata.yC; r.simdata.yT]
  ygrid = cdd.make_ygrid(y, 100)
  plot(size=plotsize)
  cdd.plot_post_density!(r.chain, ygrid)
  cdd.plot_simtruth(r.simdata.mmC, r.simdata.mmT, ygrid, label=nothing)
  savefig(joinpath(imgdir, "post-density.pdf"))
  closeall()

  # Plot gammas (zero part)
  begin
    pzero = Pzero(NC=r.simdata.NC, NT=r.simdata.NT,
                  QC=r.simdata.QC, QT=r.simdata.QT)
    pzero_inference = cdd.infer_Pzero1(pzero, 10000)
    gammaC_post = pzero_inference.distC
    gammaT_post = pzero_inference.distT
    gamma_grid = cdd.make_gamma_grid(gammaC_post, gammaT_post, a=.001)

    plot(size=plotsize)
    plot!(gamma_grid, pdf.(gammaC_post, gamma_grid), label=nothing, color=:blue)
    plot!(gamma_grid, pdf.(gammaT_post, gamma_grid), label=nothing, color=:red)
    cdd.plot_gamma_uq!(gammaC_post, color=:blue, Q=pzero.QC, N=pzero.NC,
                       truth=r.simdata.gammaC)
    cdd.plot_gamma_uq!(gammaT_post, color=:red, Q=pzero.QT, N=pzero.NT,
                       truth=r.simdata.gammaT)
    savefig(joinpath(imgdir, "gamma.pdf"))
    closeall()
  end

  # Plot F̃ᵢ
  begin 
    plot(size=plotsize)
    cdd.plot_Fi_tilde_cdf!(r.model, r.chain, pzero)
    ylims!(0, 1)
    xlabel!(L"\tilde{y}")
    ylabel!("CDF")
    savefig(joinpath(imgdir, "Fi-tilde-postmean.pdf"))
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
