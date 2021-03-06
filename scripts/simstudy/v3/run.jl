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
@time res = pmap(run, sims, on_error=identity)

# Print results status.
foreach(z -> println("$(z[1]) => $(z[2])"), zip(sims, res))

# Post process
pp_res = pmap(postprocess, sims)
# postprocess(Dict(:beta=>1, :K=>5, :snum=>1))

for snum in [1,2,3,4]
  imdir = mkpath(joinpath(Info.resultsdir_simstudy, simname, "img"))
  plot(size=plotsize)
  xlabel!("K")
  ylabel!("DIC")
  for stm in skewtmix
    dics = Float64[]
    for K in Ks
      info_path = joinpath(Info.resultsdir_simstudy, simname,
                           savename(Dict(:snum => snum, :K => K, :skewtmix => stm)),
                           "info.txt")
      model_info = open(f->read(f, String), info_path)
      _dic = parse_dic(model_info)
      append!(dics, _dic)
    end
    modelname = stm ? "Skew-t mixture" : "Normal mixture"
    plot!(Ks, dics, marker=:square, ms=8, label=modelname, legend=:topright)
  end
  savefig(joinpath(imdir, "dic_snum=$(snum).pdf"))
  closeall()
end


# Send results to S3
Util.s3sync(from="$(Info.resultsdir_simstudy)/$(simname)",
            to="$(Info.awsbucket_simstudy)/$(simname)",
            tags=`--exclude '*.nfs'`)
