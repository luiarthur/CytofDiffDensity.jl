"""
inflated: Inflated component.
rest: The rest of the distribution.
prob: Probability for the inflated component.

Note that `inflated` should not be in the support of `rest`.
"""
struct HurdleModel{A<:Real, B<:UnivariateDistribution, C<:Real} <: ContinuousUnivariateDistribution
  inflated::A
  rest::B
  prob::C
  probrest::C
end
HurdleModel(inflated, rest, prob) = HurdleModel(inflated, rest, prob, 1 - prob)

function Distributions.pdf(d::HurdleModel, x::Real)
  return x == d.inflated ? d.prob : d.probrest * pdf(d.rest, x)
end

function Distributions.logpdf(d::HurdleModel, x::Real)
  if x == d.inflated
    return log(d.prob)
  else
    return log(d.probrest) + logpdf(d.rest, x)
  end
end

function Distributions.rand(rng::Random.AbstractRNG, d::HurdleModel)
  d.prob > rand(rng) ? d.inflated : rand(rng, d.rest)
end
