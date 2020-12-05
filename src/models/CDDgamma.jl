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
  gammaC_samples = rand(dist, nsamps)

  return (gammaC_samples=gammaC_samples, gammaT_samples=gammaC_samples,
          distC=dist)
end

"""
Model1: Proportion of zeros are different.
Return posterior samples of gammaC and gammaT.
"""
function CDDg1(NC::Int, NT::Int, QC::Int, QT::Int, nsamps::Int; a::Real=1, b::Real=1)
  Nsum = NC + NT
  Qsum = QC + QT
  
  # Posterior distributions.
  distC = Beta(a + QC, b + NC - QC)
  distT = Beta(a + QT, b + NT - QT)

  # Posterior samples of gammaC and gammaT.
  gammaC_samples = rand(distC, nsamps)
  gammaT_samples = rand(distT, nsamps)

  return (gammaC_samples=gammaC_samples, gammaT_samples=gammaC_samples,
          distC=distC, distT=distT)
end
