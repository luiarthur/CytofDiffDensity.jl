println("Compile libraries on main processor..."); flush(stdout)
include("imports.jl")
println("Finished loading libraries."); flush(stdout)

using Distributed
rmprocs(workers())
addprocs(20)
println("Compile libraries on workers ..."); flush(stdout)
@everywhere include("imports.jl")
println("Finished loading libraries."); flush(stdout)

sims = dict_list(Dict(:beta => collect(0:1),
                      :K => collect(2:6),
                      :snum => collect(1:4)))

# Run simulations in parallel.
@time res = pmap(run, sims, on_error=identity)

# Print results status.
foreach(z -> println("$(z[1]) => $(z[2])"), zip(sims, res))

# Post process
pp_res = pmap(postprocess, sims)
# postprocess(Dict(:beta=>1, :K=>5, :snum=>1))

# Compute BF
for d in dict_list(Dict(:K => collect(2:6), :snum => collect(1:4)))
  d0 = copy(d); d0[:beta] = 0
  d1 = copy(d); d1[:beta] = 1
  compute_bf(d0, d1)
end

# Send results to S3
# Util.s3sync(from="$(Info.resultsdir_simstudy)/$(simname)",
#             to="$(Info.awsbucket_simstudy)/$(simname)",
#             tags=`--exclude '*.nfs'`)
