println("Compile libraries on main processor..."); flush(stdout)
include("imports.jl")
println("Finished loading libraries."); flush(stdout)

using Distributed
rmprocs(workers())
addprocs(70)
println("Compile libraries on workers ..."); flush(stdout)
@everywhere include("imports.jl")
println("Finished loading libraries."); flush(stdout)

markers = [:CD3z, :EOMES, :Perforin, :Granzyme_A, :Siglec7, :LAG3, :CD56, :CD57]
println("markers: ", markers); flush(stdout)

Ks = collect(2:9)
println("K: ", Ks)

skewtmix = [true, false]
println("skewtmix: ", skewtmix); flush(stdout)

sims = dict_list(Dict(:K => Ks, :marker => markers, :skewtmix => skewtmix))

# Run simulations in parallel.
println("Running analyses ..."); flush(stdout)
@time res = pmap(run, sims, on_error=identity)

# Print results status.
foreach(z -> println("$(z[1]) => $(z[2])"), zip(sims, res))
flush(stdout)

# Post process
println("Post process ..."); flush(stdout)
@time pp_res = pmap(postprocess, sims)

# Combine results, plot DIC.
println("Compute DIC ..."); flush(stdout)
for marker in markers
  imdir = mkpath(joinpath(Info.resultsdir_datastudy, simname, "img"))
  plot(size=plotsize)
  xlabel!("K")
  ylabel!("DIC")
  for stm in skewtmix
    dics = Float64[]
    for K in Ks
      info_path = joinpath(Info.resultsdir_datastudy, simname,
                           savename(Dict(:marker => marker, :K => K, :skewtmix => stm)),
                           "info.txt")
      model_info = open(f->read(f, String), info_path)
      _dic = parse_dic(model_info)
      append!(dics, _dic)
    end
    modelname = stm ? "Skew-t mixture" : "Normal mixture"

    # FIXME: This was hardcoded...
    pos = marker in (:LAG3, ) ? :bottomleft : :topright

    plot!(Ks, dics, marker=:square, ms=8, label=modelname, legend=pos)
  end
  savefig(joinpath(imdir, "dic_marker=$(marker).pdf"))
  closeall()
end


# Send results to S3
println("Send results to s3 ..."); flush(stdout)
Util.s3sync(from="$(Info.resultsdir_datastudy)/$(simname)",
            to="$(Info.awsbucket_datastudy)/$(simname)",
            tags=`--exclude '*.nfs'`)
