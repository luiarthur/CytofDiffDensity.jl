@testset "CDDG" begin
  Random.seed!(0)

  # NOTE: With Ni of 10000, should have around 10 it/s without loglike computation.
  # NOTE: With Ni of 1000, should have around > 100 it/s without loglike computation.
  yC = rand(SkewT(1, .6, 9, -5), 1000)
  yT = rand(SkewT(3, 1, 9, -10), 1000)

  for beta in (nothing, 0, 1)
    model = beta === nothing ? CDDG(yC, yT, 5, 1) : CDDG(yC, yT, 5, beta)
    init = MCMC.make_init_state(model)
    spl = make_sampler(model, init)

    nburn = beta === nothing ? 2 : 50
    nsamps = beta === nothing ? 2 : 100
    thin = 4
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
    cdd.printsummary(chain, metrics)

    # NOTE: This is just to check that the samples are updated every iteration.
    # Note that as discrete parameters (`lambda_i`) are present, the number of
    # unique samples may be < `nsamps`.
    if beta !== nothing
      syms = keys(init)
      for s in syms
        # println(s)
        samps = getindex.(chain, s)
        if s != :nu
          @test length(unique(samps)) > nsamps * 0.6
        else
          @test length(unique(samps)) > 2
        end
      end
    end
  end
end
