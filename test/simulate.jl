@testset "simulate" begin
  @testset "simulate_ftilde" begin
    simdata = simulate_ftilde(NC=1000, NT=1000, etaC=[.3, .3, .4], etaT=[.2,.2,.6],
                              loc=[-1,0,1], scale=[.3,.3,.3], df=[10,10,10],
                              skew=[-10,-10,-10], gammaC=0.1, gammaT=0.2)
    @test true
  end
end
