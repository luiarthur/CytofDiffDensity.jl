@testset "CDDG" begin
  # NOTE: With Ni of 10000, should have around 10 it/s without loglike computation.
  yC = rand(SkewT(1, .6, 9, -5), 10000)
  yT = rand(SkewT(3, 1, 9, -10), 10000)

  for beta in (nothing, 0, 1)
    model = beta === nothing ? CDDG(yC, yT, 5, 1) : CDDG(yC, yT, 5, beta)
    init = MCMC.make_init_state(model)
    spl = make_sampler(model, init)

    nburn = beta === nothing ? 2 : 10
    nsamps = beta === nothing ? 2 : 20
    thin = 2
    function callback(chain, state, sample, i, metrics, iterator)
      if i == 1
        metrics[:loglike_G] = Float64[]
      elseif i > nburn && mod(i, thin) == 0
        ll = loglike_G(model, state)
        append!(metrics[:loglike_G], ll)
        MCMC.ProgressBars.set_postfix(iterator, loglike=round(ll, digits=3))
      end
    end
    chain, metrics = mcmc(spl, nsamps, init=init, nburn=nburn, thin=thin, callback=callback)
    @test length(metrics[:loglike_G]) == nsamps
    @assert length(unique(chain)) == nsamps
  end
  @test true
end
