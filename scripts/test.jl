# TODO: Remove when done with this.

import Pkg; Pkg.activate("../")
using CytofDiffDensity; const cdd = CytofDiffDensity
using Distributions
using MCMC
using StatsPlots

# Normal Ordered Prior.
K = 3
m1 = MixtureModel([Normal(-2, .3), Normal(0, .2), Normal(2, .5)], [.3, .4, .3])
y = rand(m1, 500)
onm_prior = OrderedNormalMeanPrior(K, Normal(0, 6), truncated(Normal(0, 3), 0, Inf))
model1= MixSkewT(y, K, mu=onm_prior, tdist=false, skew=false)
spl = make_sampler(model1)
chain, metrics = mcmc(spl, 5000, nburn=100, thin=1)
@assert isapprox(mean.(m1.components), mean(getindex.(chain, :mu)), atol=0.2)
@assert isapprox(std.(m1.components), mean(getindex.(chain, :sigma)), atol=0.1)
plot(hcat(getindex.(chain, :mu)...)')



# SkewT Ordered Prior.
m2 = MixtureModel([SkewT(-5, .2, 10, -7),
                   SkewT(0, .3, 10, -7),
                   SkewT(5, .2, 10, -7)], [.3, .4, .3])
y = rand(m2, 500)
histogram(y, bins=50)
onm_prior = OrderedNormalMeanPrior(K, Normal(-3, 1), truncated(Normal(1, 3), 0, Inf))
model2 = MixSkewT(y, K, mu=onm_prior, skew=true, tdist=true)
spl = make_sampler(model2)
chain, metrics = mcmc(spl, 2000, nburn=0, thin=1)
mean(getindex.(chain, :mu))
plot(hcat(getindex.(chain, :mu)...)')
plot(hcat(getindex.(chain, :sigma)...)')
plot(hcat(getindex.(chain, :nu)...)')
plot(hcat(getindex.(chain, :phi)...)')
