println("Compile libraries on main processor..."); flush(stdout)
include("imports.jl")
println("Finished loading libraries."); flush(stdout)

using Distributed
rmprocs(workers())
addprocs(32)
println("Compile libraries on workers ..."); flush(stdout)
@everywhere include("imports.jl")
println("Finished loading libraries."); flush(stdout)

markers = [:CD3z, :EOMES, :Perforin, :Granzyme_A, :Siglec7, :LAG3, :CD56, :CD57]
println("markers: ", markers)

Ks = collect(2:9)
println("K: ", Ks)

skewtmix = [true, false]
println("skewtmix: ", skewtmix)

sims = dict_list(Dict(:K => Ks, :marker => markers, :skewtmix => skewtmix))

# Run simulations in parallel.
println("Running analyses ...")
@time res = pmap(run, sims, on_error=identity)

# Print results status.
foreach(z -> println("$(z[1]) => $(z[2])"), zip(sims, res))

# Post process
println("Post process ...")
@time pp_res = pmap(postprocess, sims, on_error=identity)

# Print number of small clusters
@time nc_res = pmap(sim -> print_number_of_small_clusters(sim, p=0.005),
                    cdd.MCMC.ProgressBar(sims), on_error=identity)

# Combine results, plot DIC.
println("Compute DIC ...")
foreach(marker -> plot_dic(marker, Ks, skewtmix), markers)
foreach(marker -> plot_dic(marker, Ks, skewtmix, calibrate=true), markers)

# Send results to S3
Util.s3sync(from="$(Info.resultsdir_datastudy)/$(simname)",
            to="$(Info.awsbucket_datastudy)/$(simname)",
            tags=`--exclude '*.nfs'`)
