@testset "Pzero" begin
  Random.seed!(0)
  nsamps = 50000

  @testset "same dist" begin
    log_bf = cdd.compute_log_bf(Pzero(NC=10000, NT=10000, QC=1000, QT=1000), nsamps)
    @test abs(log_bf) < 2
  end

  @testset "diff dist" begin
    log_bf = cdd.compute_log_bf(Pzero(NC=10000, NT=10000, QC=9000, QT=1000), nsamps)
    @test abs(log_bf) > 10
  end
end
