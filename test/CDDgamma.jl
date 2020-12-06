@testset "CDDgamma" begin
  nsamps = 1000
  m0 = CDDg(100, 200, 10, 25, 0, 1000)
  m1 = CDDg(100, 200, 10, 25, 1, 1000)
  @test true
end
