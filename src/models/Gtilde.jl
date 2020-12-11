struct Gtilde{Y<:AbstractVector{<:Real}, A<:Dirichlet, 
              B<:Union{Normal, <:OrderedNormalMeanPrior}, C<:Real,
              D<:LogNormal, E<:Normal, F<:Gamma} <: MCMC.Model
  yC::Y  # finite log expressions from control group
  yT::Y  # finite log expressions from treatment group
  K::Int  # number of mixture components
  beta::Bool  # model indicator
  eta::A
  mu::B
  a_omega::C
  nu::D
  psi::E
  tau::F
end

has_ordered_mu(m::Gtilde) = m.mu isa OrderedNormalMeanPrior
uses_skew(m::Gtilde) = (m.psi !== Normal(0, 0))
uses_tdist(m::Gtilde) = !(m.nu.μ > 10000 && m.nu.σ == 0)

function make_ordered_prior(yC::AbstractVector{<:Real}, yT::AbstractVector{<:Real})
  # TODO
  error("Not implemented!")
end

function Gtilde(yC, yT, K, beta; eta=Dirichlet(K, 1/K), mu=nothing,
              a_omega=2.5, nu=LogNormal(2, .5), psi=Normal(0, 3),
              tau=Gamma(0.5, 1), skew::Bool=true, tdist::Bool=true)
  skew || (psi = Normal(0, 0))
  tdist || (nu = LogNormal(20000, 0))

  if mu === nothing
    y = [yC; yT]
    mu = Normal(mean(y), std(y))
  end

  return Gtilde(yC, yT, K, Bool(beta), eta, mu, a_omega, nu, psi, tau)
end

function print_model_info(m::Gtilde)
  println("Model info:")
  for key in fieldnames(Gtilde)
    if key == :yC
      println("NC (finite): ", length(m.yC))
    elseif key == :yT
      println("NT (finite): ", length(m.yT))
    else
      println("$(key): ", getfield(m, key))
    end
  end
end

function MCMC.make_init_state(m::Gtilde)
  NC = length(m.yC)
  NT = length(m.yT)
  etaC = rand(m.eta)
  etaT = rand(m.eta)
  lambdaC = rand(Categorical(etaC), NC)
  lambdaT = rand(Categorical(etaT), NT)
  tau = rand(m.tau)
  mu = if has_ordered_mu(m)
    cumsum(rand(m.mu))
  else
    rand(m.mu, m.K)
  end
  omega = rand(InverseGamma(m.a_omega, tau), m.K)
  nu = rand(m.nu, m.K)
  psi = rand(m.psi, m.K)

  if uses_tdist(m)
    vC = rand.(Gamma.(nu[lambdaC] / 2, 2 ./ nu[lambdaC]))
    vT = rand.(Gamma.(nu[lambdaT] / 2, 2 ./ nu[lambdaT]))
  else
    vC = one.(m.yC)
    vT = one.(m.yT)
  end

  if uses_skew(m)
    zetaC = rand.(truncated.(Normal.(0, 1 ./ sqrt.(vC)), 0, Inf))
    zetaT = rand.(truncated.(Normal.(0, 1 ./ sqrt.(vT)), 0, Inf))
  else
    zetaC = zero.(m.yC)
    zetaT = zero.(m.yT)
  end

  sigma = @.scalefromaltskewt(sqrt(omega), psi)
  phi = @.skewfromaltskewt(sqrt(omega), psi)

  return (etaC=etaC, etaT=etaT, lambdaC=lambdaC, lambdaT=lambdaT, tau=tau,
          mu=mu, omega=omega, nu=nu, psi=psi, vC=vC, vT=vT, zetaC=zetaC,
          zetaT=zetaT, sigma=sigma, phi=phi)
end

Nsum(m::Gtilde) = length(m.yC) + length(m.yT)

function make_sampler(m::Gtilde; init=nothing)
  init === nothing && (init = MCMC.make_init_state(m))

  if uses_tdist(m)
    cond_nu = RWM(:nu, logprob_nu, mvarwm(init.nu), bijector=Bijectors.Log{1}())
    _update_vC = update_vC
    _update_vT = update_vT
  else
    cond_nu = Conditional(:nu, (m, s) -> fill(20000.0, m.K))
    _update_vC = (m, s) -> one.(m.yC)
    _update_vT = (m, s) -> one.(m.yT)
  end

  if uses_skew(m)
    _update_zetaC = update_zetaC
    _update_zetaT = update_zetaT
    _update_psi = update_psi
  else
    _update_zetaC = (m, s) -> zero.(m.yC)
    _update_zetaT = (m, s) -> zero.(m.yT)
    _update_psi = (m, s) -> zeros(m.K)
  end

  if has_ordered_mu(m)
    _update_mu = (m, s) -> let
      lambda = [s.lambdaC; s.lambdaT]
      v = [s.vC; s.vT]
      zeta = [s.zetaC; s.zetaT]
      y = [m.yC; m.yT]

      s2 = s.omega[lambda] ./ v
      g = y - s.psi[lambda] .* zeta
      update_mu(m.mu, s.mu, g, s2, lambda)
    end
  else
    _update_mu = update_mu
  end

  return Gibbs(m,
               Conditional(:etaC, update_etaC),
               Conditional(:etaT, update_etaT),
               Conditional(:lambdaC, update_lambdaC),
               Conditional(:lambdaT, update_lambdaT),
               Conditional(:vC, _update_vC),
               Conditional(:vT, _update_vT),
               Conditional(:zetaC, _update_zetaC),
               Conditional(:zetaT, _update_zetaT),
               Conditional(:mu, _update_mu),
               Conditional(:tau, update_tau),
               Conditional(:omega, update_omega),
               cond_nu,
               Conditional(:psi, _update_psi),
               Conditional((:sigma, :phi), update_sigma_phi))
end
