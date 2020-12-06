struct CDDgamma{A<:Real, B<:Real}
  NC::Int
  NT::Int
  QC::Int
  QT::Int
  beta::Bool
  a::A
  b::B
end

function CDDgamma(; NC, NT, QC, QT, beta, a=1.0, b=1.0)
  return CDDgamma(NC, NT, QC, QT, Bool(beta), a, b)
end

"""
Model0: Proportion of zeros are the same.
Return posterior samples of gammaC and gammaT.
"""
function infer_CDDgamma0(m::CDDgamma, nsamps::Int)
  Nsum = m.NC + m.NT
  Qsum = m.QC + m.QT
  
  # Posterior distribution.
  dist = Beta(m.a + Qsum, m.b + Nsum - Qsum)

  # Posterior samples of gammaC (and gammaT)
  gamma_samples = rand(dist, nsamps)

  return (gammaC_samples=gamma_samples,
          gammaT_samples=gamma_samples,
          distC=dist, distT=dist)
end

"""
Model1: Proportion of zeros are different.
Return posterior samples of gammaC and gammaT.
"""
function infer_CDDgamma1(m::CDDgamma, nsamps::Int)
  # Posterior distributions.
  distC = Beta(m.a + m.QC, m.b + m.NC - m.QC)
  distT = Beta(m.a + m.QT, m.b + m.NT - m.QT)

  # Posterior samples of gammaC and gammaT.
  gammaC_samples = rand(distC, nsamps)
  gammaT_samples = rand(distT, nsamps)

  return (gammaC_samples=gammaC_samples,
          gammaT_samples=gammaT_samples,
          distC=distC, distT=distT)
end

function infer(m::CDDgamma, nsamps::Int)
  if m.beta
    return infer_CDDgamma0(m::CDDgamma, nsamps)
  else
    return infer_CDDgamma1(m::CDDgamma, nsamps)
  end
end
