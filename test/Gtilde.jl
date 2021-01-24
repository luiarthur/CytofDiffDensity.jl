@testset "Gtilde" begin
  Random.seed!(0)

  # NOTE: With Ni of 10000, should have around 10 it/s without loglike computation.
  # NOTE: With Ni of 1000, should have around > 100 it/s without loglike computation.
  NC = 1100
  NT = 1150
  QC = 100
  QT = 150
  yC = rand(SkewT(1, .6, 9, -5), NC - QC)
  yT = rand(SkewT(3, 1, 9, -10), NT - QT)
  K = 5
  onm_prior = OrderedNormalMeanPrior(K, Normal(0, 6), truncated(Normal(1, .3), 0, Inf))

  for arg in [(skew=true, tdist=true, mu=onm_prior), (skew=false, tdist=false, mu=nothing)]
    for beta in (nothing, 0, 1)

      model = beta === nothing ?
              Gtilde(yC, yT, K, 1, skew=arg.skew, tdist=arg.tdist, mu=arg.mu) :
              Gtilde(yC, yT, K, beta, skew=arg.skew, tdist=arg.tdist, mu=arg.mu)

      init = MCMC.make_init_state(model)
      spl = make_sampler(model, init=init)

      print_model_info(model); @test true

      nburn = beta === nothing ? 2 : 50
      nsamps = beta === nothing ? 2 : 100
      thin = 4
      callback = make_callback(model, nburn=nburn, nsamps=nsamps, thin=thin)
      chain, metrics = mcmc(spl, nsamps, init=init, nburn=nburn, thin=thin, callback=callback)
      @test length(metrics[:loglike]) == nsamps
      @assert length(unique(chain)) == nsamps
      cdd.printsummary(chain, metrics)

      # Hellinger
      if beta == 1 || beta == 0
        @test let 
          pzero = Pzero(NC=NC, NT=NT, QC=QC, QT=QT)
          H = cdd.hellinger(chain, pzero, 5000)
          H_mean = mean(H)
          println("H: ", H_mean)
          0 < H_mean < 1
        end
      end

      # NOTE: This is just to check that the samples are updated every iteration.
      # Note that as discrete parameters (`lambda_i`) are present, the number of
      # unique samples may be < `nsamps`.
      if beta !== nothing
        syms = keys(init)
        for s in syms
          # println(s)
          samps = getindex.(chain, s)
          if s == :nu  # Metropolis
            if cdd.uses_tdist(model)
              @test length(unique(samps)) > 2
            else
              @test length(unique(samps)) == 1
            end
          elseif s in (:lambdaC, :lambdaT)  # discrete
            @test length(unique(samps)) > nsamps * 0.1
          elseif !cdd.uses_skew(model) && s in (:zetaC, :zetaT, :psi, :phi)
            @test length(unique(samps)) == 1
          elseif !cdd.uses_tdist(model) && s in (:vC, :vT)
            @test length(unique(samps)) == 1
          else
            @test length(unique(samps)) == nsamps
          end
        end
      end
    end
  end

  # TODO: Test that Ordered prior works.
end
