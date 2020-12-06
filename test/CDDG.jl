@testset "CDDG" begin
  yC = rand(SkewT(1, .6, 9, -5), 10000)
  yT = rand(SkewT(3, 1, 9, -10), 10000)
  model = CDDG(yC, yT, 5, 1)
  init = MCMC.make_init_state(model)
  spl = make_sampler(model, init)
  chain, metrics = mcmc(spl, 10, init=init, discard=50, thin=2)
  @test true
end
