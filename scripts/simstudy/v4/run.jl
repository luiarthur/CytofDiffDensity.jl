println("Compile libraries on main processor..."); flush(stdout)
include("imports.jl")
println("Finished loading libraries."); flush(stdout)

using Distributed
rmprocs(workers())
addprocs(32)
println("Compile libraries on workers ..."); flush(stdout)
@everywhere include("imports.jl")
println("Finished loading libraries."); flush(stdout)

Ks = collect(2:9)
snums = collect(1:4)
skewtmix = [true, false]

sims = dict_list(Dict(:K => Ks, :snum => snums, :skewtmix => skewtmix))

# Run simulations in parallel.
println("Fit models in parallel...")
@time res = pmap(run, sims, on_error=identity)

# Print results status.
foreach(z -> println("$(z[1]) => $(z[2])"), zip(sims, res))

# Post process
println("Post process in parallel...")
pp_res = pmap(postprocess, sims)
# postprocess(Dict(:beta=>1, :K=>5, :snum=>1))

# Print number of small clusters
sc_res = pmap(sim -> print_number_of_small_clusters(sim, p=0.005),
              cdd.MCMC.ProgressBar(sims))

# DIC
println("Compute DIC sequentially...")
foreach(snum -> plot_dic(snum, Ks, skewtmix), snums)
foreach(snum -> plot_dic(snum, Ks, skewtmix, calibrate=true), snums)


# Send results to S3
Util.s3sync(from="$(Info.resultsdir_simstudy)/$(simname)",
            to="$(Info.awsbucket_simstudy)/$(simname)",
            tags=`--exclude '*.nfs'`)
