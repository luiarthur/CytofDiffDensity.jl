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
