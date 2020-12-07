function loglike(m::MixSkewT, s::T) where T
  loc = s.mu
  scale = s.sigma
  skew = s.phi
  df = s.nu

  kernel = skewtlogpdf.(loc', scale', df', skew', m.y)
  ll = sum(logsumexp(kernel .+ log.(s.eta)', dims=2))

  return ll
end


function loglike_G(m::CDDG, s::T) where T
  loc = s.mu
  scale = s.sigma
  skew = s.phi
  df = s.nu

  ll = 0.0
  for i in (:C, :T)
    eta = grab(s, :eta, i)
    y = grab(m, :y, i)
    kernel = skewtlogpdf.(loc', scale', df', skew', y)
    ll += sum(logsumexp(kernel .+ log.(eta)', dims=2))
  end

  return ll
end

function loglike_gamma(m::CDDgamma, nsamps::Int)
  m0 = MCMC.BangBang.setpropert(m, :beta, false)
  m1 = MCMC.BangBang.setpropert(m, :beta, true)

  r0 = infer(m0, nsamps)
  r1 = infer(m1, nsamps)

  # TODO
end

# TODO: Put this in MCMC.jl?
# Returns Bayes factor in favor of model 1 (second arg)
function log_bayes_factor(ll0::AbstractArray{<:Real}, ll1::AbstractArray{<:Real})
  N0 = length(ll0)
  N1 = length(ll1)
  ll0_harmonic_mean = log(N0) - logsumexp(-ll0)
  ll1_harmonic_mean = log(N1) - logsumexp(-ll1)

  # Ratio of harmonic means.
  return ll1_harmonic_mean - ll0_harmonic_mean
end

# Compute posterior probability: P(beta=1 | data)
function compute_pm1_G(ll0::AbstractArray{<:Real}, ll1::AbstractArray{<:Real}, p::Real)
  log_bf = log_bayes_factor(ll0, ll1)
  log_prior_odds = logit(p)
  return logistic(log_bf + log_prior_odds)
end

# Concatenate samples
function cat_G_samples(chain0, chain1, ll0, ll1, p)
  pm1 = compute_pm1_G(ll0, ll1, p)
  @assert length(chain0) == length(chain1)
  nsamps = length(chain0)
  chain = [pm1 > rand() ? chain1[b] : chain0[b] for b in 1:nsamps]
end
