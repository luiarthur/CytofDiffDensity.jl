import Pkg; Pkg.activate("../../")
include(joinpath(@__DIR__, "../info.jl"))

using Util
import Random
using DrWatson
using MCMC
using CytofDiffDensity; const cdd = CytofDiffDensity
using BSON
using Distributions  # for loading results.

simname = "v1"
make_resultsdir(sim) = mkpath("$(Info.resultsdir_compare)/$(simname)/$(savename(sim))")

function make_callback(model; nburn, nsamps, thin)
  return function callback(chain, state, sample, i, metrics, iterator)
    flush(stdout)
    if i == 1
      metrics[:loglike] = Float64[]
    elseif i > nburn && mod(i, thin) == 0
      ll = loglike(model, state)
      append!(metrics[:loglike], ll)
      MCMC.ProgressBars.set_postfix(iterator, loglike=round(ll, digits=3))
    end
  end
end

function run(sim)  # redirect output?
  Random.seed!(0)
  true_dist = SkewT(2, 1, 7, -7)
  y = rand(true_dist, 1000)

  # Make results dir if needed.
  resultsdir = make_resultsdir(sim)

  # TODO: Plot and save data
  model = MixSkewT(y, sim[:K])
  init = MCMC.make_init_state(model)
  spl = make_sampler(model, init, sim[:skew], sim[:tdist])

  nburn = 300 # 3000
  nsamps = 200 # 2000
  thin = 2

  callback = make_callback(model, nburn=nburn, nsamps=nsamps, thin=thin)
  Util.redirect_stdout_to_file(joinpath(resultsdir, "log.txt")) do
    # Run MCMC.
    chain, metrics = mcmc(spl, nsamps, init=init, nburn=nburn, thin=thin,
                          callback=callback, exclude=[:v, :zeta, :lambda]);

    # Save chain and metrics.
    BSON.bson("$(resultsdir)/results.bson",
              Dict(:chain => chain, :metrics => metrics, :true_dist => true_dist,
                   :model => model, :y => y))
  end
end
