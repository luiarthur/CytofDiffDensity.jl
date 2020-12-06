@testset "SkewT" begin
  # Test it compiles.
  st = SkewT(randn()*3, rand(), rand(LogNormal(3.5, .5)), rand()*3)

  skew = -3
  scale = 2
 
  altskew = toaltskew(scale, skew)
  altscale = toaltscale(scale, skew)

  @test isapprox(scalefromaltskewt(altscale, altskew), scale)
  @test isapprox(skewfromaltskewt(altscale, altskew), skew)

  st_samples = rand(st, 100000)
end
