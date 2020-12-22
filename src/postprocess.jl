make_ygrid(y, ngrid=100) = range(minimum(y), maximum(y), length=ngrid)

function postmean_cdf(chain::AbstractVector, ygrid)
  pp = posterior_predictive(chain)
  B = length(pp.C)
  N = length(ygrid)
  cdf_C = mean(map(p -> rand(p, N) .< ygrid, pp.C))
  cdf_T = mean(map(p -> rand(p, N) .< ygrid, pp.T))
  return (C=cdf_C, T=cdf_T)
end

function posterior_predictive(chain::AbstractVector)
  function pp(c::NamedTuple, i::Symbol)
    comps = SkewT.(c.mu, c.sigma, c.nu, c.phi) 
    MixtureModel(comps, grab(c, :eta, i))
  end
  return (C=pp.(chain, :C), T=pp.(chain, :T))
end

function compute_post_density(chain::AbstractVector, ygrid)
  function compute_pdf(c::NamedTuple, i::Symbol)
    comps = SkewT.(c.mu, c.sigma, c.nu, c.phi) 
    return pdf.(MixtureModel(comps, grab(c, :eta, i)), ygrid)
  end
  return (C=compute_pdf.(chain, :C), T=compute_pdf.(chain, :T))
end


function printsummary(chain, metrics; digits=3)
  if :loglike in keys(metrics)
    println("mean loglike: ", mean(metrics[:loglike]))
    println("DIC: ", MCMC.dic(metrics[:loglike]))
  end

  sanitize(s) = round(mean(s), digits=digits)
  sanitizevec(v) = round.(mean(v), digits=digits)

  function print_uni_param(sym)
    param = getindex.(chain, sym)
    println("mean $(sym): ", sanitize(param))
  end

  function print_mv_param(sym)
    param = getindex.(chain, sym)
    println("mean $(sym): ", sanitizevec(param))
  end

  foreach(print_mv_param, [:etaC, :etaT, :mu, :sigma, :nu, :phi])
  foreach(print_uni_param, [:tau])
end


# TODO: test
function plot_post_density!(chain, ygrid; plot_mean=true, alpha=0.3)
  pdfC, pdfT = compute_post_density(chain, ygrid)
  pdfC = hcat(pdfC...)
  pdfT = hcat(pdfT...)

  pdfC_mean = mean(pdfC, dims=2)
  pdfC_lower = vec(quantiles(pdfC, 0.025, dims=2))
  pdfC_upper = vec(quantiles(pdfC, 0.975, dims=2))

  pdfT_mean = mean(pdfT, dims=2)
  pdfT_lower = vec(quantiles(pdfT, 0.025, dims=2))
  pdfT_upper = vec(quantiles(pdfT, 0.975, dims=2))

  p = let
    plot!(ygrid, pdfC_lower, fillrange=pdfC_upper, alpha=alpha, color=:blue,
         label=nothing)
    plot_mean && plot!(ygrid, pdfC_mean, color=:blue, label=nothing)
    plot!(ygrid, pdfT_lower, fillrange=pdfT_upper, alpha=alpha, color=:red,
          label=nothing)
    plot_mean && plot!(ygrid, pdfT_mean, color=:red, label=nothing)
    ylabel!("density")
  end

  return p
end


function plot_simtruth(dC, dT, ygrid; kwargs...)
  plot!(ygrid, pdf.(dC, ygrid), color=:blue, ls=:dot; kwargs...)
  plot!(ygrid, pdf.(dT, ygrid), color=:red, ls=:dot; kwargs...)
  ylabel!("density")
end


# TODO: add methods for computing and plotting CDF of F_i
function compute_Fi_tilde_cdf(gtilde::Gtilde, chain::AbstractVector,
                              pzero::Pzero; gridsize::Int=100)
  gammaC_mean = mean(infer_Pzero1(pzero, 10000).distC)
  gammaT_mean = mean(infer_Pzero1(pzero, 10000).distT)

  ygrid = make_ygrid([gtilde.yC; gtilde.yT], gridsize)
  cdfs = postmean_cdf(chain, ygrid)

  cdfC = cdfs.C * (1 - gammaC_mean) .+ gammaC_mean
  cdfT = cdfs.T * (1 - gammaT_mean) .+ gammaT_mean

  # TODO: Compute area between curves.
  area = sum(abs.(cdfC - cdfT) * (ygrid[2] - ygrid[1]))

  return (C=cdfC, T=cdfT, ygrid=ygrid, area=area)
end
function plot_Fi_tilde_cdf!(gtilde::Gtilde, chain::AbstractVector,
                            pzero::Pzero; gridsize::Int=100)
  cdfs = compute_Fi_tilde_cdf(gtilde, chain, pzero, gridsize=gridsize)
  plot!(cdfs.ygrid, cdfs.C, lw=2, color=:blue, label=nothing)
  plot!(cdfs.ygrid, cdfs.T, lw=2, color=:red, label=nothing)

  println("Area between curves: $(cdfs.area)")
end

function compute_Fi_tilde_cdf_truth(simdata::NamedTuple) end
function plot_Fi_tilde_cdf_truth!(simdata::NamedTuple) end
