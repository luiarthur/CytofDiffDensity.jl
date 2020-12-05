module CytofDiffDensity

using MCMC
using Distributions
using StatsFuns
using Random

export CDDG1, CDDG0, CDDg0, CDDg1

abstract type CDDG <: MCMC.Model end

include("models/CDDG1.jl")
include("models/CDDG0.jl")
include("models/CDDgamma.jl")

end # module
