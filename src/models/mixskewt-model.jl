struct MixSkewT{T<:AbstractVector{<:Real}, A<:Dirichlet, B<:Normal, C<:Real,
                D<:LogNormal, E<:Normal, F<:Gamma} <: MCMC.Model
  y::T
  K::Int
  eta::A
  mu::B
  a_omega::C
  nu::D
  psi::E
  tau::F
end

function MixSkewT(y, K; eta=Dirichlet(K, 1/K), mu=Normal(), a_omega=2.5,
                  nu=LogNormal(2, .5), psi=Normal(0, 3), tau=Gamma(0.5, 1))
  return MixSkewT(y, K, eta, mu, a_omega, nu, psi, tau)
end

function MCMC.make_init_state(m::MixSkewT)
  N = length(m.y)
  eta = rand(m.eta)
  lambda = rand(Categorical(eta), N)
  tau = rand(m.tau)
  mu = rand(m.mu, m.K)
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

function make_sampler(m::MixSkewT, init)
  return Gibbs(m,
               Conditional(:eta, update_eta),
               Conditional(:lambda, update_lambda),
               Conditional(:v, update_v),
               Conditional(:zeta, update_zeta),
               Conditional(:mu, update_mu),
               Conditional(:tau, update_tau),
               Conditional(:omega, update_omega),
               RWM(:nu, logprob_nu, mvarwm(init.nu), bijector=Bijectors.Log{1}()),
               Conditional(:psi, update_psi),
               Conditional((:sigma, :phi), update_sigma_phi))
end
