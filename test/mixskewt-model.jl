@testset "MixSkewT" begin
  Random.seed!(0)
  println("Test MixSkewT")

  y = rand(SkewT(1, .6, 9, -5), 1000)

  model = MixSkewT(y, 5)
  init = MCMC.make_init_state(model)
  spl = make_sampler(model, init=init)

  nburn = 50
  nsamps = 100
  thin = 2

  chain, metrics = mcmc(spl, nsamps, init=init, nburn=nburn, thin=thin)
  @test length(unique(chain)) == nsamps

  # TODO:
  # - [ ] Test that the parameters are correct.

  syms = keys(init)
  for s in syms
    # println(s)
    samps = getindex.(chain, s)
    if s != :nu
      @test length(unique(samps)) > nsamps * 0.8
    else
      @test length(unique(samps)) > 2
    end
  end
end
