struct Pzero{A<:Real, B<:Real}
  NC::Int
  NT::Int
  QC::Int
  QT::Int
  a::A
  b::B
end

function Pzero(; NC, NT, QC, QT, a=1.0, b=1.0)
  return Pzero(NC, NT, QC, QT, a, b)
end

"""
Model0: Proportion of zeros are the same.
Return posterior samples of gammaC and gammaT.
"""
function infer_Pzero0(m::Pzero, nsamps::Int)
  Nsum = m.NC + m.NT
  Qsum = m.QC + m.QT
  
  # Posterior distribution.
  dist = Beta(m.a + Qsum, m.b + Nsum - Qsum)
  
  # Posterior samples of gammaC (and gammaT)
  gamma_samples = rand(dist, nsamps)

  # Loglikelihood
  loglike = Qsum .* log.(gamma_samples) + (Nsum - Qsum) .* log1p.(-gamma_samples)

  return (gammaC_samples=gamma_samples,
          gammaT_samples=gamma_samples,
          distC=dist, distT=dist, loglike=loglike)
end

"""
Model1: Proportion of zeros are different.
Return posterior samples of gammaC and gammaT.
"""
function infer_Pzero1(m::Pzero, nsamps::Int)
  # Posterior distributions.
  distC = Beta(m.a + m.QC, m.b + m.NC - m.QC)
  distT = Beta(m.a + m.QT, m.b + m.NT - m.QT)

  # Posterior samples of gammaC and gammaT.
  gammaC_samples = rand(distC, nsamps)
  gammaT_samples = rand(distT, nsamps)

  # Loglikelihood
  loglike = let
    llC = m.QC .* log.(gammaC_samples) + (m.NC - m.QC) .* log1p.(-gammaC_samples)
    llT = m.QT .* log.(gammaT_samples) + (m.NT - m.QT) .* log1p.(-gammaT_samples)
    llC + llT
  end

  return (gammaC_samples=gammaC_samples,
          gammaT_samples=gammaT_samples,
          distC=distC, distT=distT, loglike=loglike)
end

function compute_log_bf(m::Pzero, nsamps::Int)
  m0 = infer_Pzero0(m::Pzero, nsamps)
  m1 = infer_Pzero1(m::Pzero, nsamps)

  return MCMC.log_bayes_factor(m0.loglike, m1.loglike)
end

function plot_gamma_uq!(gamma_mean::Real, gamma_lower::Real, gamma_upper::Real;
                        color, Q, truth=nothing, alpha=.6)
  if truth === nothing
    # Plot empirical mean
    vline!(Q, color=color, label=nothing, ls=:dash)
  else
    # Plot simulation truth
    vline!(truth, color=color, label=nothing, ls=:dash)
  end

  # Plot posterior mean
  vline!(gamma_mean, , color=color, alpha=alpha, label=nothing)

  # 95% CI
  vline!([gamma_lower, gamma_upper], color=color, label=nothing, ls=:dot, lw=2)

  xlabel!(L"\gamma_i")
  ylabel!("density")
end

function plot_gamma_uq!(gamma::Distribution; ci_level::Real=0.05, color, Q,
                        truth=nothing, alpha=.6)
  gamma_mean = mean(gamma)
  gamma_lower = quantile(gamma, ci_level/2)
  gamma_upper = quantile(gamma, 1 - ci_level/2)

  plot_gamma_uq!(gamma_mean, gamma_lower, gamma_upper, color=color, Q=Q,
                 truth=truth, alpha=alpha)
end

function plot_gamma_uq!(gamma::Union{Distribution, AbstractVector{<:Real}};
                        color, Q, ci_level::Real=0.05, truth=nothing, alpha=.6)
  gamma_mean = mean(gamma)
  gamma_lower = quantile(gamma, ci_level/2)
  gamma_upper = quantile(gamma, 1 - ci_level/2)
  plot_gamma_uq!(gamma_mean, gamma_lower, gamma_upper, color=color, Q=Q,
                 truth=truth, alpha=alpha)
end


function make_gamma_grid(gammaC_post::Distribution, gammaT_post::Distribution;
                         a=0.01, grid_length=100)
  gamma_lower = min(quantile(gammaC_post, a), quantile(gammaT_post, a))
  gamma_upper = max(quantile(gammaC_post, 1-a), quantile(gammaT_post, 1-a))
  return range(gamma_lower, gamma_upper, length=grid_length)
end
