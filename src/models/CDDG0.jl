struct CDDG0{Y<:AbstractVector{<:Real}, A<:Dirichlet, B<:Normal, C<:Real,
             D<:LogNormal, E<:Normal, F<:Gamma} <: MCMC.Model
  yC::Y  # finite log expressions from control group
  yT::Y  # finite log expressions from treatment group
  K::Int  # number of mixture components
  eta::A
  mu::B
  a_omega::C
  nu::D
  psi::E
  tau::F
end
function CDDG0(yC, yT, K; eta=Dirichlet(K, 1/K), mu=Normal(0, 3), a_omega=2.5,
               nu=LogNormal(2, .5), psi=Normal(0, 3), tau=Gamma(0.5, 1))
  return CDDG0(yC, yT, K, eta, mu, a_omega, nu, psi, tau)
end

function MCMC.make_init_state(m::CDDG0)
  # TODO: Return a NamedTuple
end

function update_mu(m::CDDG0, s::T) where T
  # TODO: Implement updates.
  return s.mu
end

# TODO: Implement Gibbs sampler.
