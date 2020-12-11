println("Test OrderedNormalMeanPrior")

@testset "OrderedNormalMeanPrior" begin
  @testset "GMM (ordered prior)" begin
    Random.seed!(0)
    K = 3
    mm = MixtureModel([Normal(-2, .3), Normal(0, .2), Normal(2, .5)], [.3, .4, .3])
    y = rand(mm, 500)

    onm_prior = OrderedNormalMeanPrior(K, Normal(0, 6), truncated(Normal(1, .3), 0, Inf))
    model = MixSkewT(y, K, mu=onm_prior, tdist=false, skew=false)
    print_model_info(model)
    spl = make_sampler(model)

    nsamps = 1000
    chain, metrics = mcmc(spl, nsamps, nburn=2000, thin=1)
    @test length(unique(chain)) == nsamps
    @test isapprox(mean.(mm.components), mean(getindex.(chain, :mu)), atol=0.2)
    @test isapprox(std.(mm.components), mean(getindex.(chain, :sigma)), atol=0.1)
  end

  @testset "MixSkewT (ordered prior)" begin
    Random.seed!(0)
    K = 3
    mm = MixtureModel([SkewT(-5, .2, 10, -7),
                       SkewT(0, .3, 10, -7),
                       SkewT(5, .2, 10, -7)], [.3, .4, .3])
    y = rand(mm, 1000)

    # NOTE:
    # - Prior for psi needs to be loose (e.g. Normal(0, 3). Normal(-2, 1) is
    #   too strict and skew may be extremely negative in posterior.)
    # - OrderedNormalMeanPrior needs to somewhat loose: if too strict, then
    #   difficult to escape bad local modes.
    onm_prior = OrderedNormalMeanPrior(K, Normal(0, 3), truncated(Normal(1, 3), 0, Inf))
    model = MixSkewT(y, K, mu=onm_prior, psi=Normal(0, 3), eta=Dirichlet(K, 1))
    print_model_info(model)

    # Find best seed.
    seed, min_loglike = find_good_seed(model, nsamps=100, nburn=200, thin=2, seeds=41:60)
    println("Best seed: $seed ($min_loglike)")

    # Define callback.
    nsamps, nburn, thin = (2000, 3000, 2)
    callback = make_callback(model, nsamps=nsamps, nburn=nburn, thin=thin)

    # Run official model.
    Random.seed!(seed)
    spl = make_sampler(model)
    chain, metrics = mcmc(spl, nsamps, nburn=nburn, thin=thin, callback=callback)

    # Tests.
    @test length(unique(chain)) == nsamps
    @test isapprox(getfield.(mm.components, :loc), mean(getindex.(chain, :mu)), atol=0.1)
    @test isapprox(getfield.(mm.components, :scale), mean(getindex.(chain, :sigma)), atol=0.1)
    @test let
      # Test that the 95% CI contains truth.
      nu_post = hcat(getindex.(chain, :nu)...)
      nu_true = getfield.(mm.components, :df)
      nu_lower = vec(MCMC.quantiles(nu_post, .025, dims=2))
      nu_upper = vec(MCMC.quantiles(nu_post, .975, dims=2))
      all(nu_lower .< nu_true .< nu_upper)
    end
    @test isapprox(getfield.(mm.components, :skew), mean(getindex.(chain, :phi)), atol=4)
    @test isapprox(mm.prior.p, mean(getindex.(chain, :eta)), atol=.1)

    # TODO: check that the 95% CI contains truth.
  end
end
