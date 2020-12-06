using CytofDiffDensity
using Test

using Distributions
using StatsFuns
import Random

@testset "MCMC" begin
  include("skewt.jl")
  include("CDDgamma.jl")
end
