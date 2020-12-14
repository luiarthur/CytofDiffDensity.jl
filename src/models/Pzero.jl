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
