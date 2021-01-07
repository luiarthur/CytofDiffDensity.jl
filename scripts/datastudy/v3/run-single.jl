println("Compile libraries on main processor..."); flush(stdout)
include("imports.jl")
println("Finished loading libraries."); flush(stdout)

println("ARGS: ", ARGS)

marker = Symbol(ARGS[1])
println("marker: ", marker); flush(stdout)

K = parse(Int, ARGS[2])
println("K: ", K)

stm = parse(Bool, ARGS[3])
println("skewtmix: ", stm); flush(stdout)

sims = dict_list(Dict(:K => K, :marker => marker, :skewtmix => stm))

# Run simulations in parallel.
println("Running analyses ..."); flush(stdout)
@time run(sims[1])

# Post process
println("Post process ..."); flush(stdout)
@time postprocess(sims[1])

# Print number of small clusters
@time print_number_of_small_clusters(sims[1], p=0.02)
