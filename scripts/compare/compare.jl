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

# Post process.
using StatsPlots

# Plot settings
ENV["GKSwstype"] = "nul"  # For StatsPlots to plot in background only.
plotsize = (450, 450)
Plots.scalefontsizes()
Plots.scalefontsizes(1.5)

# Get all DICs.
dics = Dict(
  map(sim -> savename(sim) => 
      dic(BSON.load(joinpath(make_resultsdir(sim), "results.bson"))[:metrics][:loglike]),
      sims))

sims1 = filter(sim -> occursin("skew=true_tdist=true", savename(sim)), sims)
sims2 = filter(sim -> occursin("skew=true_tdist=false", savename(sim)), sims)
sims3 = filter(sim -> occursin("skew=false_tdist=true", savename(sim)), sims)
sims4 = filter(sim -> occursin("skew=false_tdist=false", savename(sim)), sims)

imgdir = mkpath(joinpath(Info.resultsdir_compare, simname, "img"))
plot(1:5, map(sim -> getindex(dics, savename(sim)), sims4), palette=:tab10,
     lw=5, marker=:, ms=8, label="n-mix")
plot!(1:5, map(sim -> getindex(dics, savename(sim)), sims3), palette=:tab10,
      lw=5, marker=:, ms=8, label="t-mix")
plot!(1:5, map(sim -> getindex(dics, savename(sim)), sims2), palette=:tab10,
      lw=5, marker=:, ms=8, label="skew-n mix")
plot!(1:5, map(sim -> getindex(dics, savename(sim)), sims1), palette=:tab10,
      lw=5, marker=:, ms=8, label="skew-t mix")
plot!(size=plotsize)
xlabel!("K")
ylabel!("DIC")
savefig(joinpath(imgdir, "dic-compare.pdf"))
closeall()
