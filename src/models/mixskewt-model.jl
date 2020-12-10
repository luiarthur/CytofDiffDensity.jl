struct MixSkewT{T<:AbstractVector{<:Real}, A<:Dirichlet,
                B<:Union{Normal, <:OrderedNormalMeanPrior},
                C<:Real, D<:LogNormal, E<:Normal, F<:Gamma} <: MCMC.Model
  y::T
  K::Int
  eta::A
  mu::B
  a_omega::C
  nu::D
  psi::E
  tau::F
end

has_ordered_mu(m::MixSkewT) = m.mu isa OrderedNormalMeanPrior

function MixSkewT(y, K; eta=Dirichlet(K, 1/K), mu=Normal(), a_omega=2.5,
                  nu=LogNormal(2, .5), psi=Normal(0, 3), tau=Gamma(0.5, 1))
  return MixSkewT(y, K, eta, mu, a_omega, nu, psi, tau)
end

function MCMC.make_init_state(m::MixSkewT)
  N = length(m.y)
  eta = rand(m.eta)
  lambda = rand(Categorical(eta), N)
  tau = rand(m.tau)
  mu = if has_ordered_mu(m)
    cumsum(rand(m.mu))
  else
    rand(m.mu, m.K)
  end
  omega = rand(InverseGamma(m.a_omega, tau), m.K)
  nu = rand(m.nu, m.K)
  psi = rand(m.psi, m.K)
  v = rand.(Gamma.(nu[lambda] / 2, 2 ./ nu[lambda]))
  zeta = rand.(truncated.(Normal.(0, 1 ./ sqrt.(v)), 0, Inf))

  sigma = @.scalefromaltskewt(sqrt(omega), psi)
  phi = @.skewfromaltskewt(sqrt(omega), psi)

  return (eta=eta, lambda=lambda, tau=tau, mu=mu, omega=omega, nu=nu, psi=psi,
          v=v, zeta=zeta, sigma=sigma, phi=phi)
end

# TODO: Remove the following if no problems arise.
#
# function _make_sampler(m::MixSkewT, init)
#   return Gibbs(m,
#                Conditional(:lambda, update_lambda),
#                Conditional(:v, update_v),
#                Conditional(:zeta, update_zeta),
#                Conditional(:eta, update_eta),
#                Conditional(:mu, update_mu),
#                Conditional(:tau, update_tau),
#                Conditional(:omega, update_omega),
#                Conditional(:psi, update_psi),
#                RWM(:nu, logprob_nu, mvarwm(init.nu), bijector=Bijectors.Log{1}()),
#                Conditional((:sigma, :phi), update_sigma_phi))
# end

"""
    make_sampler(m::MixSkewT; init, skew::Bool=true, tdist::Bool=true)

Note that when `skew=false` and `tdist=false`, then we have a Gaussian mixture
model.
"""
function make_sampler(m::MixSkewT; init=nothing, skew::Bool=true, tdist::Bool=true)
  init === nothing && (init = MCMC.make_init_state(m))
  if tdist
    cond_nu = RWM(:nu, logprob_nu, mvarwm(init.nu), bijector=Bijectors.Log{1}())
    _update_v = update_v
  else
    cond_nu = Conditional(:nu, (m, s) -> fill(20000.0, m.K))
    _update_v = (m, s) -> one.(m.y)
  end

  if skew
    _update_zeta = update_zeta
    _update_psi = update_psi
  else
    _update_zeta = (m, s) -> zero.(m.y)
    _update_psi = (m, s) -> zeros(m.K)
  end

  if has_ordered_mu(m)
    _update_mu = (m, s) -> update_mu(m.mu, s.mu, m.y, s.omega[s.lambda], s.lambda)
  else
    _update_mu = update_mu
  end

  return Gibbs(m,
               Conditional(:lambda, update_lambda),
               Conditional(:v, _update_v),
               Conditional(:zeta, _update_zeta),
               Conditional(:eta, update_eta),
               Conditional(:mu, _update_mu),
               Conditional(:tau, update_tau),
               Conditional(:omega, update_omega),
               Conditional(:psi, _update_psi),
               cond_nu,
               Conditional((:sigma, :phi), update_sigma_phi))
end


function print_model_info(m::MixSkewT)
  println("Model info:")
  for key in fieldnames(MixSkewT)
    if key == :y
      println("N: ", length(y))
    else
      println("$(key): ", getfield(m, key))
    end
  end
end
