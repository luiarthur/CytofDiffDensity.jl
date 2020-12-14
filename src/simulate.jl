function simulate_gtilde(; NC, NT, etaC, etaT, loc, scale, df, skew)
  K = length(loc)

  @assert K == length(etaC) == length(etaT)
  @assert K == length(scale) == length(df) == length(skew)

  mmC = MixtureModel(SkewT.(loc, scale, df, skew), etaC)
  mmT = MixtureModel(SkewT.(loc, scale, df, skew), etaT)

  yC = rand(mmC, NC)
  yT = rand(mmT, NT)

  return (yC=yC, yT=yT, etaC=etaC, etaT=etaT, loc=loc, scale=scale, df=df,
          skew=skew, mmC=mmC, mmT=mmT)
end


function simulate_ftilde(; NC, NT, etaC, etaT, loc, scale, df, skew, gammaC, gammaT)
  QC = sum(rand(NC) .> gammaC)
  QT = sum(rand(NT) .> gammaT)

  NC_finite = NC - QC
  NT_finite = NT - QT

  pzero = (gammaC=gammaC, gammaT=gammaT, QC=QC, QT=QT, NC=NC, NT=NT)

  gtilde = simulate_gtilde(NC=NC_finite, NT=NT_finite, etaC=etaC, etaT=etaT, loc=loc,
                           scale=scale, df=df, skew=skew)

  return merge(pzero, gtilde)
end
