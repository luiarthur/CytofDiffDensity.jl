using CytofDiffDensity
using Test

using MCMC
using Distributions
using StatsFuns
import Random

@testset "MCMC" begin
  include("skewt.jl")
  include("CDDgamma.jl")
  include("CDDG.jl")
end
