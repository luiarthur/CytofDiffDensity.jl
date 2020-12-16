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

# TODO: Combine results
for snum in [1,2,3,4]
  imdir = mkpath(joinpath(Info.resultsdir_simstudy, simname, "img"))
  dics = Float64[]
  for K in MCMC.ProgressBar([2,3,4,5,6])
    cr = combine_results(snum, K)
    append!(dics, cr.dic)
  end
  open(joinpath(imdir, "dics_snum=$(snum).txt"), "w") do io
    write(io, "DICs: $(dics) \n")
  end
  # Model comparisons
  plot(2:6, dics, marker=:square, ms=8, label=nothing)
  xlabel!("K")
  ylabel!("DIC")
  savefig(joinpath(imdir, "dic_snum=$(snum).pdf"))
  closeall()
end


# Send results to S3
Util.s3sync(from="$(Info.resultsdir_simstudy)/$(simname)",
            to="$(Info.awsbucket_simstudy)/$(simname)",
            tags=`--exclude '*.nfs'`)
