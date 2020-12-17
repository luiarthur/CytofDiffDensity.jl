# TODO:
# Modify:
# - [ ] Make K more interesting across scenarios (but `K <= 5`).
# - [ ] make Ni = 100_000
# - [ ] Change gamma in the different scenarios.

function scenarios(n::Int; seed=nothing, Ni=1000)
  seed == nothing || Random.seed!(seed)

  if n == 0  # TEST CASE!
    return scenarios(1, seed=seed, Ni=100)
  elseif n == 1
    return simulate_ftilde(NC=Ni, NT=Ni,
                           etaC=[.25, .75, 0, 0],
                           etaT=[.1, .1, .8, 0], 
                           loc = [-1.5, 3.5, 5.1, 5],
                           scale = [1.6, 1.76, 1.76, 1.6],
                           df = [12, 10, 10, 15.],
                           skew = [0, -10, -10, 0.],
                           gammaC=0.1, gammaT=0.1)
  elseif n == 2
    return simulate_ftilde(NC=Ni, NT=Ni,
                           etaC=[.1, .1, .5, .3], 
                           etaT=[.1, .1, .8, 0], 
                           loc = [-1.5, 3.5, 5.1, 4.3],
                           scale = [1.6, 1.76, 1.76, 1.6],
                           df = [12, 10, 10, 15.],
                           skew = [0, -10, -10, -11.],
                           gammaC=0.1, gammaT=0.1)
  elseif n == 3
    return simulate_ftilde(NC=Ni, NT=Ni,
                           etaC=[.05, .05, .5, .4], 
                           etaT=[.05, .05, .5, .4], 
                           loc = [-1.5, 3.5, 5.1, 4.3],
                           scale = [1.6, 1.76, 1.76, 1.6],
                           df = [12, 10, 10, 15.],
                           skew = [0, -10, -10, -11.],
                           gammaC=0.1, gammaT=0.2)
  elseif n == 4
    return simulate_ftilde(NC=Ni, NT=Ni,
                           etaC=[.05, .05, .5, .4], 
                           etaT=[.05, .05, .5, .4], 
                           loc = [-1.5, 3.5, 5.1, 4.3],
                           scale = [1.6, 1.76, 1.76, 1.6],
                           df = [12, 10, 10, 15.],
                           skew = [0, -10, -10, -11.],
                           gammaC=0.15, gammaT=0.15)
  elseif n == 5  # a test case
    return simulate_ftilde(NC=Ni, NT=Ni,
                           etaC=[0, 1., 0, 0],
                           etaT=[0, 0, 1., 0],
                           loc = [-1.5, 3.5, 5.1, 4.3],
                           scale = [1.6, 1.76, 1.76, 1.6],
                           df = [12, 10, 10, 15.],
                           skew = [0, -10, -10, -11.],
                           gammaC=0.1, gammaT=0.1)
  else
    throw(ArgumentError("n=$(n) was provided. But n must be an integer between 0 and 5."))
  end
end
