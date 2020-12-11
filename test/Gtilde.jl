@testset "Gtilde" begin
  Random.seed!(0)

  # NOTE: With Ni of 10000, should have around 10 it/s without loglike computation.
  # NOTE: With Ni of 1000, should have around > 100 it/s without loglike computation.
  yC = rand(SkewT(1, .6, 9, -5), 1000)
  yT = rand(SkewT(3, 1, 9, -10), 1000)

  for beta in (nothing, 0, 1)
    model = beta === nothing ? Gtilde(yC, yT, 5, 1) : Gtilde(yC, yT, 5, beta)
    init = MCMC.make_init_state(model)
    spl = make_sampler(model, init=init)

    print_model_info(model); @test true

    nburn = beta === nothing ? 2 : 50
    nsamps = beta === nothing ? 2 : 100
    thin = 4
    callback = make_callback(model, nburn=nburn, nsamps=nsamps, thin=thin)
    chain, metrics = mcmc(spl, nsamps, init=init, nburn=nburn, thin=thin, callback=callback)
    @test length(metrics[:loglike]) == nsamps
    @assert length(unique(chain)) == nsamps
    cdd.printsummary(chain, metrics)

    # NOTE: This is just to check that the samples are updated every iteration.
    # Note that as discrete parameters (`lambda_i`) are present, the number of
    # unique samples may be < `nsamps`.
    if beta !== nothing
      syms = keys(init)
      for s in syms
        # println(s)
        samps = getindex.(chain, s)
        if s != :nu
          @test length(unique(samps)) > nsamps * 0.6
        else
          @test length(unique(samps)) > 2
        end
      end
    end
  end

  # TODO: Test that Ordered prior works.
end
