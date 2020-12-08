println("Compile on main node ...")
include("imports.jl")
using Distributed

rmprocs(workers())
addprocs(20)  # 20
println("Compile on workers ...")
@everywhere include("imports.jl")

# Simulations to run.
sims = dict_list(Dict(:skew => [true, false], :tdist => [true, false], :K => collect(1:5)))

# Run simulations in parallel.
@time res = pmap(run, sims, on_error=identity)

# Send results to S3
Util.s3sync(from="$(Info.resultsdir_compare)/$(simname)",
            to="$(Info.awsbucket_compare)/$(simname)",
            tags=`--exclude '*.nfs'`)

# Postproessing
include("postprocess.jl")
