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
Plots.scalefontsizes(1.4)

simname = "v2"
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
  mu_prior = cdd.make_ordered_prior(y, sim.K)  # ordered prior for mu.
  # mu_prior = Normal(mean(y), std(y))  # independent priors for mu_k.

  if sim.skewtmix
    # SkewT mixture
    model = Gtilde(data.yC, data.yT, sim.K, true, mu=mu_prior, tau=Gamma(1, 1),
                   skew=true, tdist=true)
  else
    # Normal mixture
    model = Gtilde(data.yC, data.yT, sim.K, true, mu=mu_prior, tau=Gamma(1, 1),
                   skew=false, tdist=false)
  end
  exclude = [:vC, :vT, :zetaC, :zetaT, :lambdaC, :lambdaT]
                 
  init = MCMC.make_init_state(model)
  spl = make_sampler(model, init=init)

  nburn, nsamps, thin = (40000, 4000, 2)
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
  println("Data info:")
  foreach(s -> println(s, ": ", getindex(r.data, s)), [:QC, :QT, :NC, :NT])
  println()

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

    # TODO: Print hellinger distance
    H = cdd.hellinger(r.chain, pzero, 10000, verbose=0)
    H_mean = mean(H)
    H_lower = quantile(H, .025)
    H_upper = quantile(H, .975)
    println("H mean: ", H_mean)
    println("H lower: ", H_lower)
    println("H upper: ", H_upper)

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
    cdd.plot_Fi_tilde_cdf!(r.model, r.chain, pzero, N=2000)
    ylims!(0, 1)
    xlabel!(L"\tilde{y}")
    ylabel!("CDF")
    savefig(joinpath(imgdir, "Fi-tilde-postmean.pdf"))
    closeall()

    # Plot Fi
    plot(size=plotsize)
    cdd.plot_Fi_tilde_cdf!(r.model, r.chain, pzero, exponentiate=true, N=2000)
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


function print_number_of_small_clusters(sim; p=0.02)
  # Make results dir if needed.
  resultsdir = make_resultsdir(sim)

  Util.redirect_stdout_to_file(joinpath(resultsdir, "small_clusters.txt")) do
    _print_number_of_small_clusters(sim, resultsdir, p=p)
  end
end


function _print_number_of_small_clusters(sim, resultsdir; p=0.02)
  r = (; BSON.load(joinpath(resultsdir, "results.bson"))...)
  num_small_clusters = cdd.count_small_clusters(r.chain, p=p)
  println(num_small_clusters)
end


function plot_dic(marker, Ks, skewtmix; calibrate=false)
  imdir = mkpath(joinpath(Info.resultsdir_datastudy, simname, "img"))

  function init_plot(; half=false, use_xlabel=true)
    if half
      plot(size=(plotsize[1], plotsize[2] / 2))
    else
      plot(size=plotsize)
    end

    if calibrate
      use_xlabel && xlabel!("number of small components")
    else
      use_xlabel && xlabel!("K")
    end

    ylabel!("DIC")
  end

  calibrate || init_plot()
  calibrate_plots = []
  for stm in skewtmix
    dics = Float64[]
    num_small_clusters = Float64[]
    for K in Ks
      expname = savename(Dict(:marker => marker, :K => K, :skewtmix => stm))
      resdir = joinpath(Info.resultsdir_datastudy, simname, expname)

      # DIC
      info_path = joinpath(resdir, "info.txt")
      model_info = open(f->read(f, String), info_path)
      _dic = parse_dic(model_info)
      append!(dics, _dic)

      # Calibration
      if calibrate
        num_small_clusters_path = joinpath(resdir, "small_clusters.txt")
        _num_small_clusters = parse(Float64, open(f->read(f, String), num_small_clusters_path))
        append!(num_small_clusters, _num_small_clusters)
      end
    end
    modelname = stm ? "Skew-t mixture" : "Normal mixture"


    if calibrate
      init_plot(use_xlabel=!stm)
      _p = plot!(num_small_clusters, dics, label=modelname,
                 legend=:topright, lw=3, alpha=.7)
      append!(calibrate_plots, [_p])
      annotate!(num_small_clusters, dics, Ks)
      savefig(joinpath(imdir, "calibrate_marker=$(marker)_skewtmix=$(stm).pdf"))
      closeall()
    else
      plot!(Ks, dics, marker=:square, ms=8, label=modelname, legend=:topright)
    end
  end

  if calibrate 
    plot(calibrate_plots[1], calibrate_plots[2], layout=@layout [a; b])
    savefig(joinpath(imdir, "calibrate_merged_marker=$(marker).pdf"))
    closeall()
  end

  calibrate || savefig(joinpath(imdir, "dic_marker=$(marker).pdf"))
  calibrate || closeall()
end

# TODO:
# - [ ] Implement `resume`.
