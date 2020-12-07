println("Compile on main node ...")
include("imports.jl")
using Distributed

rmprocs(workers())
addprocs(2)  # 20
println("Compile on workers ...")
@everywhere include("imports.jl")

# Simulations to run.
sims = dict_list(Dict(:skew => [true, false], :tdist => [true, false], :K => collect(1:5)))

# Run simulations in parallel.
@time res = pmap(run, sims, on_error=identity)

# Send results to S3
cdd.s3ync(from="$(Info.resultsdir_compare)/$(simname)",
          to="$(Info.awsbucket_compare)/$(simname)",
          tags=`--exclude '*.nfs'`)

# Post process.
# using StatsPlots
# result = BSON.load(joinpath(make_resultsdir(sims[1]), "results.bson"))
# plot(result[:metrics][:loglike])
# 
# dics = Dict(
#   map(sim -> savename(sim) => 
#       dic(BSON.load(joinpath(make_resultsdir(sim), "results.bson"))[:metrics][:loglike]),
#       sims))
# 
# sims1 = filter(sim -> occursin("skew=true_tdist=true", savename(sim)), sims)
# sims2 = filter(sim -> occursin("skew=true_tdist=false", savename(sim)), sims)
# sims3 = filter(sim -> occursin("skew=false_tdist=true", savename(sim)), sims)
# sims4 = filter(sim -> occursin("skew=false_tdist=false", savename(sim)), sims)
# 
# plot(map(sim -> getindex(dics, savename(sim)), sims1))
# plot(map(sim -> getindex(dics, savename(sim)), sims2))
# plot(map(sim -> getindex(dics, savename(sim)), sims3))
# plot(map(sim -> getindex(dics, savename(sim)), sims4))
