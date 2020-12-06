"""
Model0: Proportion of zeros are the same.
Return posterior samples of gammaC and gammaT.
"""
function CDDg0(NC::Int, NT::Int, QC::Int, QT::Int, nsamps::Int; a::Real=1, b::Real=1)
  Nsum = NC + NT
  Qsum = QC + QT
  
  # Posterior distribution.
  dist = Beta(a + Qsum, b + Nsum - Qsum)

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
function CDDg1(NC::Int, NT::Int, QC::Int, QT::Int, nsamps::Int; a::Real=1, b::Real=1)
  # Posterior distributions.
  distC = Beta(a + QC, b + NC - QC)
  distT = Beta(a + QT, b + NT - QT)

  # Posterior samples of gammaC and gammaT.
  gammaC_samples = rand(distC, nsamps)
  gammaT_samples = rand(distT, nsamps)

  return (gammaC_samples=gammaC_samples,
          gammaT_samples=gammaT_samples,
          distC=distC, distT=distT)
end


function CDDg(NC, NT, QC, QT, beta, nsamps; a=1, b=1)
  if Bool(beta)
    return CDDg1(NC, NT, QC, QT, nsamps, a=a, b=b)
  else
    return CDDg0(NC, NT, QC, QT, nsamps, a=a, b=b)
  end
end
