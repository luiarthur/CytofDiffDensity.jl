@testset "HurdleModel" begin
  Random.seed!(0)
  m1 = HurdleModel(-Inf, MixtureModel(SkewT.([0, 3], 2, 4, -10)), .3)
  m2 = HurdleModel(-Inf, MixtureModel(Normal.([0, 3])), .3)
  for m in (m1, m2)
    rand(m)
    @test mean([rand(m) for _ in 1:10000] .=== -Inf) ≈ m1.prob atol=0.01
    @test pdf(m, 10.0) ≈ exp(logpdf(m, 10)) atol=0.01
    @test pdf(m, -Inf) ≈ exp(logpdf(m, -Inf)) atol=0.01
  end
end
