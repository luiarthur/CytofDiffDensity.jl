"""
Computes hellinger distance between two distributions `f` and `g`.  Evaluates
densities `f` and `g` over a grid `xs` evenly spaced between `min` and `max`,
with length of `length`.
"""
function hellinger(f, g, min, max; length=1000)
  xs = range(min, max, length=length)
  d = xs[2] - xs[1]
  kernel = sum(d * sqrt.(f.(xs) .* g.(xs)))
  return sqrt(1 - kernel)
end

"""
Computes hellinger distance between two distributions `f` and `g`.  Evaluates
log densities `log_f` and `log_g` over a grid `xs` evenly spaced between `min`
and `max`, with length of `length`. This might be more numerically stable than
`hellinger`.
"""
function hellinger_logpdf(log_f, log_g, min, max; length=1000)
  xs = range(min, max, length=length)
  log_d = log(xs[2] - xs[1])
  log_kernel = logsumexp(log_d .+ (log_f.(xs) + log_g.(xs)) / 2)
  return sqrt(1 - exp(log_kernel))
end

function loglike(m::MixSkewT, s::T) where T
  loc = s.mu
  scale = s.sigma
  skew = s.phi
  df = s.nu

  kernel = skewtlogpdf.(loc', scale', df', skew', m.y)
  ll = sum(logsumexp(kernel .+ log.(s.eta)', dims=2))

  return ll
end


function loglike(m::Gtilde, s::T) where T
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

function loglike(m::Pzero, nsamps::Int)
  m0 = MCMC.BangBang.setpropert(m, :beta, false)
  m1 = MCMC.BangBang.setpropert(m, :beta, true)

  r0 = infer(m0, nsamps)
  r1 = infer(m1, nsamps)

  # TODO
end

# TODO: Test.
# Compute posterior probability: P(beta=1 | data)
function compute_pm1(::Gtilde, ll0::AbstractArray{<:Real}, ll1::AbstractArray{<:Real}, p::Real)
  log_bf = MCMC.log_bayes_factor(ll0, ll1)
  log_prior_odds = logit(p)
  return logistic(log_bf + log_prior_odds)
end

# TODO: Test.
# Concatenate samples
function cat_samples(m::Gtilde, chain0, chain1, ll0, ll1, p)
  pm1 = compute_pm1(m, ll0, ll1, p)
  @assert length(chain0) == length(chain1)
  nsamps = length(chain0)
  chain = [pm1 > rand() ? chain1[b] : chain0[b] for b in 1:nsamps]
end
