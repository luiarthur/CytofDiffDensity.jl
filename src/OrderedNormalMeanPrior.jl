abstract type ConjugatePrior end

struct OrderedNormalMeanPrior{T} <: ConjugatePrior
  priors::T
  function OrderedNormalMeanPrior(priors::T) where T
    # Assert that the first element is valid.
    priors[1] isa Normal || 
    priors[1] isa Truncated{<:Normal} || 
    error("Invalid prior. First element must be a Normal or Truncated{<:Normal}!")

    # Assert that all other elements are Truncated{<:Normal}
    if length(priors) > 1
      for k in 2:length(priors)
        priors[k] isa Truncated{<:Normal} ||
        error("Invalid prior. Element $(k) must be a Truncated{<:Normal}!")
      end
    end

    return new{T}(priors)
  end
end

function Distributions.rand(rng::Random.AbstractRNG, d::OrderedNormalMeanPrior)
  return rand.(rng, d.priors)
end

function OrderedNormalMeanPrior(K::Int, comp1::Union{Normal, Truncated{<:Normal}},
                                compk::Truncated{<:Normal})
  K >= 2 || error("`compk` was provided, but provided `K` was < 2!")
  return OrderedNormalMeanPrior([comp1; fill(compk, K - 1)])
end

isreal(p::OrderedNormalMeanPrior) = p.priors[1] isa Normal
ncomponents(p::OrderedNormalMeanPrior) = length(p.priors)

get_untruncated_mean(d::Normal) = d.μ
get_untruncated_mean(d::Truncated{<:Normal}) = get_untruncated_mean(d.untruncated)
get_untruncated_std(d::Normal) = d.σ
get_untruncated_std(d::Truncated{<:Normal}) = get_untruncated_std(d.untruncated)

function get_untruncated_mean_std(d::Union{Normal,Truncated{<:Normal}})
  return (get_untruncated_mean(d), get_untruncated_std(d))
end
   
# TODO: Test
function update(p::OrderedNormalMeanPrior, curr::AbstractVector{<:Real},
                y::AbstractVector{<:Real}, sigmasq::AbstractVector{<:Real},
                lambda::AbstractVector{<:Integer}) 
  K = ncomponents(p)

  # New value.
  next = copy(curr)

  # Update `next` sequentially.
  for k in 1:K
    # Set of indices such that lambda[i] > k
    Ck = findall(lambda .>= k)

    # Prior (mean, sd) of untruncated Normal
    mk, sk = get_untruncated_mean_std(p.priors[k])

    # Compute new variance for untruncated Normal
    sigmasq_k = sigmasq[Ck]
    vnew = 1 / (1/sk^2 +  sum(1 ./ sigmasq_k))

    # Compute new mean for untruncated Normal
    next[k] = 0
    mu = cumsum(next)
    lam = lambda[Ck]
    mnew = vnew * (mk / sk^2 + sum((y[Ck] - mu[lam]) ./ sigmasq_k))

    d_untruncated = Normal(mnew, sqrt(vnew))

    if p.priors[k] isa Normal
      d = d_untruncated
    else
      d = truncated(d_untruncated, p.priors[k].lower, p.priors[k].upper)
    end

    next[k] = rand(d)
  end

  return next
end

to_mu(delta) = cumsum(delta)
to_delta(mu) = [mu[1]; diff(mu)]

function update_mu(p::OrderedNormalMeanPrior, curr::AbstractVector{<:Real},
                   y::AbstractVector{<:Real}, sigmasq::AbstractVector{<:Real},
                   lambda::AbstractVector{<:Integer})
  delta = to_delta(curr)
  new_delta = update(p, delta, y, sigmasq, lambda)
  return to_mu(new_delta)
end
