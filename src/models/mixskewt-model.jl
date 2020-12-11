struct MixSkewT{T<:AbstractVector{<:Real}, A<:Dirichlet,
                B<:Union{Normal, <:OrderedNormalMeanPrior},
                C<:Real, D<:LogNormal, E<:Normal, F<:Gamma} <: MCMC.Model
  y::T
  K::Int
  eta::A
  mu::B
  a_omega::C
  nu::D  # If want GMM, set `nu = LogNormal(m, 0.0)`, where `m > 10000`.
  psi::E  # If want no skew, set `psi = Normal(0, 0)`.
  tau::F
end

has_ordered_mu(m::MixSkewT) = m.mu isa OrderedNormalMeanPrior

function MixSkewT(y, K; eta=Dirichlet(K, 1/K), mu=Normal(), a_omega=2.5,
                  nu=LogNormal(2, .5), psi=Normal(0, 3), tau=Gamma(0.5, 1),
                  skew::Bool=true, tdist::Bool=true)
  skew || (psi = Normal(0, 0))
  tdist || (nu = LogNormal(20000, 0))
  return MixSkewT(y, K, eta, mu, a_omega, nu, psi, tau)
end

uses_skew(m::MixSkewT) = (m.psi !== Normal(0, 0))
uses_tdist(m::MixSkewT) = !(m.nu.μ > 10000 && m.nu.σ == 0)

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

  if uses_tdist(m)
    v = rand.(Gamma.(nu[lambda] / 2, 2 ./ nu[lambda]))
  else
    v = one.(m.y)
  end

  if uses_skew(m)
    zeta = rand.(truncated.(Normal.(0, 1 ./ sqrt.(v)), 0, Inf))
  else
    zeta = zero.(m.y)
  end

  sigma = @.scalefromaltskewt(sqrt(omega), psi)
  phi = @.skewfromaltskewt(sqrt(omega), psi)

  return (eta=eta, lambda=lambda, tau=tau, mu=mu, omega=omega, nu=nu, psi=psi,
          v=v, zeta=zeta, sigma=sigma, phi=phi)
end

"""
    make_sampler(m::MixSkewT; init, skew::Bool=true, tdist::Bool=true)

Note that when `skew=false` and `tdist=false`, then we have a Gaussian mixture
model.
"""
function make_sampler(m::MixSkewT; init=nothing)
  init === nothing && (init = MCMC.make_init_state(m))

  if uses_tdist(m)
    cond_nu = RWM(:nu, logprob_nu, mvarwm(init.nu), bijector=Bijectors.Log{1}())
    _update_v = update_v
  else
    cond_nu = Conditional(:nu, (m, s) -> fill(20000.0, m.K))
    _update_v = (m, s) -> one.(m.y)
  end

  if uses_skew(m)
    _update_zeta = update_zeta
    _update_psi = update_psi
  else
    _update_zeta = (m, s) -> zero.(m.y)
    _update_psi = (m, s) -> zeros(m.K)
  end

  if has_ordered_mu(m)
    _update_mu = (m, s) -> let
      s2 = s.omega[s.lambda] ./ s.v
      g = m.y - s.psi[s.lambda] .* s.zeta
      update_mu(m.mu, s.mu, g, s2, s.lambda)
    end
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
      println("N: ", length(m.y))
    else
      println("$(key): ", getfield(m, key))
    end
  end
end


function make_callback(model::MCMC.Model; nburn, nsamps, thin)
  return function callback(chain, state, sample, i, metrics, iterator)
    flush(stdout)
    if i == 1
      metrics[:loglike] = Float64[]
    elseif i > nburn && mod(i, thin) == 0
      ll = loglike(model, state)
      append!(metrics[:loglike], ll)
      if iterator isa MCMC.ProgressBar
        MCMC.ProgressBars.set_postfix(iterator, loglike=round(ll, digits=3))
      end
    end
  end
end


"""
Finding a good starting value is important when using OrderedNormalMeanPrior.
"""
function find_good_seed(model::MCMC.Model; seeds=1:10, nsamps=1, nburn=200, thin=1)
  # Preallocate results in Dict{seed => mean(loglike)}.
  results = Dict((seed, 0.0) for seed in seeds)
  callback = make_callback(model, nsamps=nsamps, nburn=nburn, thin=thin)

  for seed in MCMC.ProgressBar(seeds)
    Random.seed!(seed)
    spl = make_sampler(model)
    chain, metrics = mcmc(spl, nsamps, nburn=nburn, thin=thin,
                          callback=callback, progress=false)
    results[seed] = mean(metrics[:loglike])
  end

  best_seed = argmax(results)
  return best_seed, results[best_seed]
end
