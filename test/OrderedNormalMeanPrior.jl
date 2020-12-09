println("Test OrderedNormalMeanPrior")

@testset "OrderedNormalMeanPrior" begin
  @testset "update" begin
    K = 3
    delta_prior = OrderedNormalMeanPrior(K, Normal(0, 6), truncated(Normal(1, .3), 0, Inf))
    mm = MixtureModel([Normal(-2, .3), Normal(0, .2), Normal(2, .5)], [.3, .4, .3])
    y = rand(mm, 500)

    function update_mu(m::MixSkewT, s::T) where T
      return cdd.update_mu(delta_prior, s.mu, m.y, s.omega[s.lambda], s.lambda)
    end

    model = MixSkewT(y, K)
    spl = Gibbs(model,
                Conditional(:lambda, cdd.update_lambda),
                Conditional(:v, (m, s) -> one.(m.y)),
                Conditional(:zeta, (m, s) -> zero.(m.y)),
                Conditional(:eta, cdd.update_eta),
                Conditional(:mu, update_mu),
                Conditional(:tau, cdd.update_tau),
                Conditional(:omega, cdd.update_omega),
                Conditional(:psi, (m, s) -> zeros(m.K)),
                Conditional(:nu, (m, s) -> fill(20000.0, m.K)),
                Conditional((:sigma, :phi), cdd.update_sigma_phi))

    nsamps = 1000
    chain, metrics = mcmc(spl, nsamps, nburn=2000, thin=1)
    @test length(unique(chain)) == nsamps
    @test isapprox(mean.(mm.components), mean(getindex.(chain, :mu)), atol=1e-1)
    @test isapprox(std.(mm.components), mean(getindex.(chain, :sigma)), atol=1e-1)
  end
end
