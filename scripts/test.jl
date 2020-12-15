# TODO: Remove when done with this.

import Pkg; Pkg.activate("../")
using CytofDiffDensity; const cdd = CytofDiffDensity
using Distributions
using MCMC
using StatsPlots

function make_ordered_prior(y, K; s=nothing, upper=.9, lower=.1)
  s === nothing && (s = std(y) / K)
  comp1_mean = quantile(y, lower)
  compk_mean = (quantile(y, upper) - quantile(y, lower)) / (K - 1)
  return OrderedNormalMeanPrior(K,
                                Normal(comp1_mean, s),
                                truncated(Normal(compk_mean, s), 0, Inf))
end

# Normal Ordered Prior.
K = 3
m1 = MixtureModel([Normal(-2, .3), Normal(0, .2), Normal(1, .4)], [.3, .4, .3])
y = rand(m1, 1000)
histogram(y, label=nothing, bins=100)
onm_prior = make_ordered_prior(y, K)
model1= MixSkewT(y, K, mu=onm_prior, tdist=false, skew=false)
spl = make_sampler(model1)
chain, metrics = mcmc(spl, 5000, nburn=100, thin=1)
@assert isapprox(mean.(m1.components), mean(getindex.(chain, :mu)), atol=0.2)
@assert isapprox(std.(m1.components), mean(getindex.(chain, :sigma)), atol=0.1)
plot(hcat(getindex.(chain, :mu)...)')



# SkewT Ordered Prior.
K = 3
m2 = MixtureModel([SkewT(-5, 1, 10, -7),
                   SkewT(0, 2, 10, -12),
                   SkewT(3, 3, 10, -10)], [.3, .4, .3])
y = rand(m2, 1000)

# m2 = MixtureModel([SkewT(-5, .2, 10, -7),
#                    SkewT(0, .3, 10, -7),
#                    SkewT(5, .2, 10, -7)], [.3, .4, .3])
# y = rand(m2, 1000)

histogram(y, bins=100, label=nothing)
onm_prior = make_ordered_prior(y, K, s=1)
model2 = MixSkewT(y, K, mu=onm_prior, skew=true, tdist=true)
spl = make_sampler(model2)
chain, metrics = mcmc(spl, 1000, nburn=0, thin=3)
mean(getindex.(chain, :mu))
mean(getindex.(chain, :eta))
plot(hcat(getindex.(chain, :mu)...)', label=nothing)
hline!(getfield.(m2.components, :loc), label=nothing, ls=:dot, color=:black)
plot(hcat(getindex.(chain, :sigma)...)')
plot(hcat(getindex.(chain, :nu)...)')
plot(hcat(getindex.(chain, :phi)...)')
plot(hcat(getindex.(chain, :eta)...)')

