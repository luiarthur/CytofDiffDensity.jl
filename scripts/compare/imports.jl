ENV["GKSwstype"] = "nul"  # For StatsPlots to plot in background only.

import Pkg; Pkg.activate("../../")
include(joinpath(@__DIR__, "../info.jl"))

using Util
import Random
using DrWatson
using MCMC
using CytofDiffDensity; const cdd = CytofDiffDensity
using BSON
using Distributions  # for loading results.
using StatsPlots

# Plot settings
plotsize = (450, 450)
Plots.scalefontsizes()
Plots.scalefontsizes(1.5)

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
  spl = make_sampler(model, init, skew=sim[:skew], tdist=sim[:tdist])

  nburn = 3000
  nsamps = 2000
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


"""
Plots histogram of data and true data density (dotted line).
"""
function plot_simdata(y, true_dist; grid_length=100)
  ygrid = range(minimum(y), maximum(y), length=grid_length)
  histogram(y, label=nothing, normalize=true, la=0)
  plot!(ygrid, pdf.(true_dist, ygrid), label=nothing, color=:black, ls=:dot)
  ylabel!("density")
end


"""
Computes posterior density.
"""
function compute_post_density(chain, ygrid)
  function compute_pdf(c::NamedTuple)
    comps = SkewT.(c.mu, c.sigma, c.nu, c.phi) 
    return pdf.(MixtureModel(comps, c.eta), ygrid)
  end
  return compute_pdf.(chain)
end

"""
Plot posterior density and simtruth.
"""
function plot_simtruth_with_post(true_dist, chain, ygrid; alpha=0.3)
  plot(ygrid, pdf.(true_dist, ygrid), label=nothing, color=:black)
  ylabel!("density")

  post_dens = hcat(compute_post_density(chain, ygrid)...)
  pdf_lower = vec(MCMC.quantiles(post_dens, 0.025, dims=2))
  pdf_upper = vec(MCMC.quantiles(post_dens, 0.975, dims=2))

  plot!(ygrid, pdf_lower, fillrange=pdf_upper, alpha=alpha, color=:blue,
        label=nothing)
end

function postprocess(sim)
  r = let
    d = BSON.load(joinpath(make_resultsdir(sim), "results.bson"))
    (; d...)  # Dict -> NamedTuple
  end

  imgdir = mkpath(joinpath(Info.resultsdir_compare, simname, savename(sim), "img"))

  ygrid = cdd.make_ygrid(r.y)
  plot_simtruth_with_post(r.true_dist, r.chain, ygrid)
  plot!(size=plotsize)
  savefig(joinpath(imgdir, "post-density.pdf"))
  closeall()

  function plot_uniparam(sym)
    param = getindex.(r.chain, sym)
    p1 = plot(param, label=nothing)
    p2 = histogram(param, label=nothing, normalize=true)
    plot(p1, p2, layout=(1, 2), size=(1000, 500))
    savefig(joinpath(imgdir, "$(sym).pdf"))
    closeall()
  end

  function plot_mvparam(sym)
    param = hcat(getindex.(r.chain, sym)...)'
    size(param, 2) == 1 && (param = vec(param))
    p1 = plot(param, label=nothing)
    p2 = boxplot(param, label=nothing)
    plot(p1, p2, layout=(1, 2), size=(1000, 500))
    savefig(joinpath(imgdir, "$(sym).pdf"))
    closeall()
  end

  # Plot trace plots and boxplots for multivariate parameters.
  foreach(plot_mvparam, [:mu, :sigma, :nu, :phi, :eta])

  # Plot trace plots and boxplots for univariate parameters.
  foreach(plot_uniparam, [:tau])
end
