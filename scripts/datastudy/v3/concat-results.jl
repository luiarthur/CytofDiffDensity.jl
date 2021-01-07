println("Compile libraries on main processor..."); flush(stdout)
include("imports.jl")
println("Finished loading libraries."); flush(stdout)

using Distributed
rmprocs(workers())
addprocs(20)
println("Compile libraries on workers ..."); flush(stdout)
@everywhere include("imports.jl")
println("Finished loading libraries."); flush(stdout)

markers = [:CD3z, :EOMES, :Perforin, :Granzyme_A, :Siglec7, :LAG3, :CD56, :CD57]
println("markers: ", markers); flush(stdout)

Ks = collect(2:9)
println("K: ", Ks)

skewtmix = [true, false]
println("skewtmix: ", skewtmix); flush(stdout)

# Jobs 
sims = dict_list(Dict(:K => Ks, :marker => markers, :skewtmix => skewtmix))

# Combine results, plot DIC.
println("Compute DIC ..."); flush(stdout)
foreach(marker -> plot_dic(marker, Ks, skewtmix), markers)
foreach(marker -> plot_dic(marker, Ks, skewtmix, calibrate=true), markers)

# Send results to S3
println("Send results to s3 ..."); flush(stdout)
Util.s3sync(from="$(Info.resultsdir_datastudy)/$(simname)",
            to="$(Info.awsbucket_datastudy)/$(simname)",
            tags=`--exclude '*.nfs'`)
