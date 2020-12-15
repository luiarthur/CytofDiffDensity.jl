# TODO: Remove when done with this.

import Pkg; Pkg.activate("../")
using CytofDiffDensity; const cdd = CytofDiffDensity
using Distributions
using MCMC
using StatsPlots
import Random

function compute_post_density(chain::AbstractVector, ygrid)
  function compute_pdf(c::NamedTuple)
    comps = SkewT.(c.mu, c.sigma, c.nu, c.phi) 
    return pdf.(MixtureModel(comps, c.eta), ygrid)
  end
  return compute_pdf.(chain)
end


# Normal Ordered Prior.
K = 3
m1 = MixtureModel([Normal(-2, .3), Normal(0, .2), Normal(1, .4)], [.3, .4, .3])
y = rand(m1, 1000)
histogram(y, label=nothing, bins=100)
onm_prior = cdd.make_ordered_prior(y, K)
model1= MixSkewT(y, K, mu=onm_prior, tdist=false, skew=false)
spl = make_sampler(model1)
chain, metrics = mcmc(spl, 5000, nburn=100, thin=1)
@assert isapprox(mean.(m1.components), mean(getindex.(chain, :mu)), atol=0.2)
@assert isapprox(std.(m1.components), mean(getindex.(chain, :sigma)), atol=0.1)
plot(hcat(getindex.(chain, :mu)...)')


# SkewT Ordered Prior.
K = 3
Random.seed!(1)
m2 = MixtureModel([SkewT(-5, 1, 10, -7),
                   SkewT(0, 2, 10, -12),
                   SkewT(3, 3, 10, -10)], [.3, .4, .3])
y = rand(m2, 10000)
histogram(y, bins=100, label=nothing)

Random.seed!(0)
# onm_prior = cdd.make_ordered_prior(y, K, s=(maximum(y)-minimum(y))/2K)
onm_prior = cdd.make_ordered_prior(y, K)
model2 = MixSkewT(y, K, mu=onm_prior, skew=true, tdist=true)
spl = make_sampler(model2)
nburn, nsamps, thin = (3000, 1000, 1)
callback = make_callback(model2, nburn=nburn, nsamps=nsamps, thin=thin)
chain, metrics = mcmc(spl, nsamps, nburn=nburn, thin=thin, callback=callback)
plot(metrics[:loglike])

mean(getindex.(chain, :mu))
mean(getindex.(chain, :eta))
plot(hcat(getindex.(chain, :mu)...)', label=nothing)
hline!(getfield.(m2.components, :loc), label=nothing, ls=:dot, color=:black)
plot(hcat(getindex.(chain, :sigma)...)')
hline!(getfield.(m2.components, :scale), label=nothing, ls=:dot, color=:black)
plot(hcat(getindex.(chain, :nu)...)')
hline!(getfield.(m2.components, :df), label=nothing, ls=:dot, color=:black)
plot(hcat(getindex.(chain, :phi)...)')
hline!(getfield.(m2.components, :skew), label=nothing, ls=:dot, color=:black)
plot(hcat(getindex.(chain, :eta)...)')
hline!(m2.prior.p, label=nothing, ls=:dot, color=:black)


begin
  ygrid = cdd.make_ygrid(y)
  dens = compute_post_density(chain, ygrid)
  plot(ygrid, pdf.(m2, ygrid), label=nothing, color=:black)
  dens_lower = vec(MCMC.quantiles(hcat(dens...), 0.025, dims=2))
  dens_upper = vec(MCMC.quantiles(hcat(dens...), 0.975, dims=2))
  plot!(ygrid, dens_lower, fillrange=dens_upper, label=nothing, alpha=0.3)
  histogram!(y, normalize=true, la=0, alpha=.3)
end

