using CytofDiffDensity; const cdd = CytofDiffDensity
using Test

using MCMC
using Distributions
using StatsFuns
import Random

@testset "MCMC" begin
  include("hurdle.jl")
  include("simulate.jl")
  include("skewt.jl")
  include("mixskewt-model.jl")
  include("Pzero.jl")
  include("Gtilde.jl")
  include("OrderedNormalMeanPrior.jl")
end
