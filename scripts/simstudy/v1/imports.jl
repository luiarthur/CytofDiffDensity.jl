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
using Distributions  # required for `BSON.load`
using StatsPlots
using StatsFuns

# Plot settings
plotsize = (450, 450)
Plots.scalefontsizes()
Plots.scalefontsizes(1.5)

simname = "v1"
make_resultsdir(sim) = mkpath("$(Info.resultsdir_simstudy)/$(simname)/$(savename(sim))")

function run(sim)
  sim = (; sim...)

  # Simulate data.
  simdata = scenarios(sim.snum, Ni=10000, seed=1)

  # Make results dir if needed.
  resultsdir = make_resultsdir(sim)

  # Create Priors.
  y = [simdata.yC; simdata.yT]

  # Plot and save data
  imgdir = mkpath(joinpath(resultsdir, "img"))
  plot(size=plotsize)
  histogram!(simdata.yC, normalize=true, label=nothing, la=0, alpha=.3, color=:blue)
  histogram!(simdata.yT, normalize=true, label=nothing, la=0, alpha=.3, color=:red)
  ygrid = cdd.make_ygrid(y, 100)
  plot!(ygrid, pdf.(simdata.mmC, ygrid), label=nothing, color=:blue)
  plot!(ygrid, pdf.(simdata.mmT, ygrid), label=nothing, color=:red)
  savefig(joinpath(imgdir, "data.pdf"))
  closeall()

  # TODO: Does this help?
  mu_prior = cdd.make_ordered_prior(y, sim.K,
                                    s=(maximum(y) - minimum(y)) / (2 * sim.K))

  Random.seed!(0)
  model = Gtilde(simdata.yC, simdata.yT, sim.K, sim.beta, mu=mu_prior)
  exclude = [:vC, :vT, :zetaC, :zetaT, :lambdaC, :lambdaT]
                 
  init = MCMC.make_init_state(model)
  spl = make_sampler(model, init=init)

  nburn, nsamps, thin = (10000, 4000, 2)
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

  # TODO: Plot posterior density for C and T.
  y = [r.simdata.yC; r.simdata.yT]
  ygrid = cdd.make_ygrid(y, 100)
  plot(size=plotsize)
  cdd.plot_post_density!(r.chain, ygrid)
  cdd.plot_simtruth(r.simdata.mmC, r.simdata.mmT, ygrid, label=nothing)
  savefig(joinpath(imgdir, "post-density.pdf"))
  closeall()
end


function postprocess(sim)
  # Make results dir if needed.
  resultsdir = make_resultsdir(sim)

  Util.redirect_stdout_to_file(joinpath(resultsdir, "info.txt")) do
    _postprocess(sim, resultsdir)
  end
end


function load_results(snum, K, beta)
  d = Dict(:snum => snum, :K => K, :beta => beta)
  resdir = make_resultsdir(d)
  return (; BSON.load(joinpath(resdir, "results.bson"))...)
end


function combine_results(snum::Integer, K::Integer; p::Real=0.1)
  (0 < p < 1) || error("0 < p < 1 required!")

  # Results dir 
  resdir = make_resultsdir(Dict(:snum => snum, :K => K))
  imgdir = mkpath(joinpath(resdir, "img"))

  # Load results
  r0 = load_results(snum, K, 0)
  r1 = load_results(snum, K, 1)
  length(r0.chain) == length(r1.chain) || error("Chain lengths need to be equal!")

  # Loglikelihoods
  ll0 = r0.metrics[:loglike]
  ll1 = r1.metrics[:loglike]
  length(ll0) == length(ll1) || error("Loglikelihood lengths need to be equal!")

  # Number of MCMC samples
  B = length(ll0)

  # Log BF for Gtilde
  log_bf_gtilde = MCMC.log_bayes_factor(ll0, ll1)

  # Log bf for PZero
  pzero = Pzero(NC=r0.simdata.NC, NT=r0.simdata.NT,
                QC=r0.simdata.QC, QT=r0.simdata.QT)
  log_bf_pzero = cdd.compute_log_bf(pzero, B)

  # Compute P(β=1 | data)
  pm1 = logistic(log_bf_gtilde + log_bf_pzero + logit(p))

  # Combine chains and loglikelihoods.
  beta = pm1 .> rand(B)
  chain = [beta[b] ? r1.chain[b] : r0.chain[b] for b in 1:B]
  loglike = [beta[b] ? ll1[b] : ll0[b] for b in 1:B]
  m0_gamma = cdd.infer_Pzero0(pzero, B)
  m1_gamma = cdd.infer_Pzero1(pzero, B)
  gammaC = [beta[b] ? m1_gamma.gammaC_samples[b] : m0_gamma.gammaC_samples[b] for b in 1:B]
  gammaT = [beta[b] ? m1_gamma.gammaT_samples[b] : m0_gamma.gammaT_samples[b] for b in 1:B]

  # Plot gamma
  plot(size=plotsize)
  histogram!(gammaC, normalize=true, alpha=0.3, color=:blue)
  histogram!(gammaT, normalize=true, alpha=0.3, color=:red)
  savefig(joinpath(imgdir, "gamma.pdf"))
  closeall()
  
  # Plot Gtilde
  y = [r0.simdata.yC; r0.simdata.yT]
  ygrid = range(quantile(y, .05), quantile(y, .95), length=100)
  plot(size=plotsize)
  cdd.plot_post_density!(chain, ygrid)
  cdd.plot_simtruth(r0.simdata.mmC, r0.simdata.mmT, ygrid, label=nothing)
  savefig(joinpath(imgdir, "post-density.pdf"))
  closeall()

  # DIC
  merged_dic = MCMC.dic(loglike)

  # Write results
  open(joinpath(resdir, "metrics.txt"), "w") do io
    write(io, "log BF(Gtilde) in favor of model 1: $(log_bf_gtilde) \n")
    write(io, "log BF(PZero) in favor of model 1: $(log_bf_pzero) \n")
    write(io, "p: $(p) \n")
    write(io, "pm1: $(pm1) \n")
    write(io, "DIC: $(merged_dic) \n")
  end

  return (chain=chain, loglike=loglike, dic=merged_dic)
end
