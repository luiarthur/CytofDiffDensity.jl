struct Gtilde{Y<:AbstractVector{<:Real}, A<:Dirichlet, B<:Normal, C<:Real,
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

function Gtilde(yC, yT, K, beta; eta=Dirichlet(K, 1/K), mu=nothing,
              a_omega=2.5, nu=LogNormal(2, .5), psi=Normal(0, 3),
              tau=Gamma(0.5, 1))
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
  mu = rand(m.mu, m.K)
  omega = rand(InverseGamma(m.a_omega, tau), m.K)
  nu = rand(m.nu, m.K)
  psi = rand(m.psi, m.K)
  vC = rand.(Gamma.(nu[lambdaC] / 2, 2 ./ nu[lambdaC]))
  vT = rand.(Gamma.(nu[lambdaT] / 2, 2 ./ nu[lambdaT]))
  zetaC = rand.(truncated.(Normal.(0, 1 ./ sqrt.(vC)), 0, Inf))
  zetaT = rand.(truncated.(Normal.(0, 1 ./ sqrt.(vT)), 0, Inf))

  sigma = @.scalefromaltskewt(sqrt(omega), psi)
  phi = @.skewfromaltskewt(sqrt(omega), psi)

  return (etaC=etaC, etaT=etaT, lambdaC=lambdaC, lambdaT=lambdaT, tau=tau,
          mu=mu, omega=omega, nu=nu, psi=psi, vC=vC, vT=vT, zetaC=zetaC,
          zetaT=zetaT, sigma=sigma, phi=phi)
end

function make_sampler(m::Gtilde, init)
  return Gibbs(m,
               Conditional(:etaC, update_etaC),
               Conditional(:etaT, update_etaT),
               Conditional(:lambdaC, update_lambdaC),
               Conditional(:lambdaT, update_lambdaT),
               Conditional(:vC, update_vC),
               Conditional(:vT, update_vT),
               Conditional(:zetaC, update_zetaC),
               Conditional(:zetaT, update_zetaT),
               Conditional(:mu, update_mu),
               Conditional(:tau, update_tau),
               Conditional(:omega, update_omega),
               RWM(:nu, logprob_nu, mvarwm(init.nu), bijector=Bijectors.Log{1}()),
               Conditional(:psi, update_psi),
               Conditional((:sigma, :phi), update_sigma_phi))
end
