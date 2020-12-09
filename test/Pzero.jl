@testset "Pzero" begin
  nsamps = 1000
  m0 = infer(Pzero(NC=100, NT=200, QC=10, QT=25, beta=0), 1000)
  m1 = infer(Pzero(NC=100, NT=200, QC=10, QT=25, beta=1), 1000)
  @test true
end