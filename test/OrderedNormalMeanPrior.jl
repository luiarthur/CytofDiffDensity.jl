println("Test OrderedNormalMeanPrior")

@testset "OrderedNormalMeanPrior" begin
  @testset "update" begin
    Random.seed!(0)
    K = 3
    mm = MixtureModel([Normal(-2, .3), Normal(0, .2), Normal(2, .5)], [.3, .4, .3])
    y = rand(mm, 500)

    onm_prior = OrderedNormalMeanPrior(K, Normal(0, 6), truncated(Normal(1, .3), 0, Inf))
    model = MixSkewT(y, K, mu=onm_prior)
    cdd.print_model_info(model)
    spl = make_sampler(model, tdist=false, skew=false)

    nsamps = 1000
    chain, metrics = mcmc(spl, nsamps, nburn=2000, thin=1)
    @test length(unique(chain)) == nsamps
    @test isapprox(mean.(mm.components), mean(getindex.(chain, :mu)), atol=0.2)
    @test isapprox(std.(mm.components), mean(getindex.(chain, :sigma)), atol=0.1)
  end
end
