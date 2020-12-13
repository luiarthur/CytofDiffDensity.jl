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
  mu_prior = let
    tn_mean = (quantile(y, .9) - mean(y)) / sim.K
    OrderedNormalMeanPrior(sim.K,
                           Normal(quantile(y, .1), 1),
                           truncated(Normal(tn_mean), 0, Inf))
  end

  Random.seed!(0)
  model = Gtilde(simdata.yC, simdata.yT, sim.K, sim.beta, mu=mu_prior)
                 
  init = MCMC.make_init_state(model)
  spl = make_sampler(model, init=init)

  nburn, nsamps, thin = (4000, 4000, 1)
  callback = make_callback(model, nburn=nburn, nsamps=nsamps, thin=thin)

  Util.redirect_stdout_to_file(joinpath(resultsdir, "log.txt")) do
    # Run MCMC.
    chain, metrics = mcmc(spl, nsamps, init=init, nburn=nburn, thin=thin,
                          callback=callback, exclude=[:v, :zeta, :lambda]);

    # Save chain and metrics.
    BSON.bson(joinpath(resultsdir, "results.bson"),
              Dict(:chain => chain, :metrics => metrics, :simdata => simdata,
                   :model => model))
  end
end
