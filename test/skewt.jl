@testset "SkewT" begin
  @testset "scalefromaltskewt, skewfromaltskewt" begin
    skew = -3
    scale = 2
    altskew = toaltskew(scale, skew)
    altscale = toaltscale(scale, skew)

    @test isapprox(scalefromaltskewt(altscale, altskew), scale)
    @test isapprox(skewfromaltskewt(altscale, altskew), skew)
  end

  @testset "Inference" begin
    Random.seed!(0)
    st = SkewT(2, 1, 7, -10)
    y = rand(st, 1000)
    model = MixSkewT(y, 1) 
    print_model_info(model)
    nsamps, nburn, thin = (2000, 2000, 2)
    callback = make_callback(model, nsamps=nsamps, nburn=nburn, thin=thin)
    spl = make_sampler(model)
    @time chain, metrics = mcmc(spl, nsamps, nburn=nburn, thin=thin, callback=callback)

    # Tests
    @test length(unique(chain)) == nsamps
    @test isapprox(st.loc, mean(getindex.(chain, :mu))[1], atol=0.1)
    @test isapprox(st.scale, mean(getindex.(chain, :sigma))[1], atol=0.1)
    @test isapprox(st.df, mean(getindex.(chain, :nu))[1], atol=2)
    @test isapprox(st.skew, mean(getindex.(chain, :phi))[1], atol=3)
  end
end
