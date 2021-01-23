function update_eta(m::MixSkewT, s::S) where S
  anew = copy(collect(m.eta.alpha))
  foreach(lam -> anew[lam] +=1, s.lambda)
  return rand(Dirichlet(anew))
end

function update_lambda(m::MixSkewT, s::S) where S
  logeta = log.(s.eta)

  lambda = [let
              loc = s.mu + s.psi * s.zeta[n]
              scale = sqrt.(s.omega) / sqrt(s.v[n])
              logmix = normlogpdf.(loc, scale, m.y[n]) + logeta
              all(s.nu .< 10000) && (logmix += gammalogpdf.(s.nu/2, 2 ./ s.nu, s.v[n]))
              MCMC.wsample_logprob(logmix)
            end for n in eachindex(m.y)]

  return lambda
end

function update_v(m::MixSkewT, s::S) where S
  nu = s.nu[s.lambda]
  omega = s.omega[s.lambda]
  psi = s.psi[s.lambda]
  mu = s.mu[s.lambda]

  shape = nu/2 .+ 1
  rate = (nu + s.zeta.^2 + ((m.y - mu - psi .* s.zeta) .^ 2) ./ omega) / 2

  return @.rand(Gamma(shape, 1 / rate))
end

function update_zeta(m::MixSkewT, s::S) where S
  mu = s.mu[s.lambda]
  omega = s.omega[s.lambda]
  psi = s.psi[s.lambda]

  vnew = @. 1 / (s.v + (psi ^ 2) * s.v / omega)
  mnew = @. vnew * s.v * psi * (m.y - mu) / omega

  return @.rand(truncated(Normal(mnew, sqrt(vnew)), 0, Inf))
end

function update_mu(m::MixSkewT, s::S) where S
  m_mu, s_mu = params(m.mu)
  vkernels = zero.(s.mu)
  mkernels = zero.(s.mu)

  for i in (:C, :T)
    for n in eachindex(s.lambda)
      k = s.lambda[n]
      vkernels[k] += s.v[n]
      mkernels[k] += (m.y[n] - s.psi[k] * s.zeta[n]) * s.v[n]
    end
  end

  vnew = 1 ./ (s_mu^-2 .+ vkernels ./ s.omega)
  mnew = vnew .* (m_mu/(s_mu^2) .+ mkernels ./ s.omega)
  
  return randn(m.K) .* sqrt.(vnew) + mnew
end

function update_tau(m::MixSkewT, s::S) where S
  a, b = params(m.tau)
  new_shape = a + m.K * m.a_omega
  new_rate = b + sum(1 ./ s.omega)
  return rand(Gamma(new_shape, 1 / new_rate))
end

function update_sigma(m::MixSkewT, s::S) where S
  return @. scalefromaltskewt(sqrt(s.omega), s.psi)
end

function update_phi(m::MixSkewT, s::S) where S
  return @. skewfromaltskewt(sqrt(s.omega), s.psi)
end

function update_sigma_phi(m::MixSkewT, s::S) where S
  return (sigma=update_sigma(m, s), phi=update_phi(m, s))
end

function update_omega(m::MixSkewT, s::S) where S
  a, b = m.a_omega, s.tau
  akernel = zero.(s.omega)
  bkernel = zero.(s.omega)

  for n in eachindex(m.y)
    k = s.lambda[n]
    akernel[k] += 1
    bkernel[k] += s.v[n] * (m.y[n] - s.mu[k] - s.psi[k] * s.zeta[n])^2
  end

  anew = a .+ akernel / 2
  bnew = b .+ bkernel / 2

  return @.rand(InverseGamma(anew, bnew))
end

function logprob_nu(m::MixSkewT, s::S, nu::AbstractVector{<:Real}) where S
  logprior = sum(logpdf.(m.nu, nu))
  loglike = let
    nu = nu[s.lambda]
    sum(gammalogpdf.(nu / 2, 2 ./ nu, s.v))
  end
  return logprior + loglike
end

function update_psi(m::MixSkewT, s::S) where S
  m_psi, s_psi = params(m.psi)
  vkernel = zero.(s.psi)
  mkernel = zero.(s.psi)

  # Update kernels
  for n in eachindex(s.v)
    k = s.lambda[n]
    vkernel[k] += s.zeta[n]^2 * s.v[n] / s.omega[k]
    mkernel[k] += s.zeta[n] * (m.y[n] - s.mu[k]) * s.v[n] / s.omega[k]
  end

  vnew = 1 ./ (s_psi^-2 .+ vkernel)
  mnew = vnew .* (m_psi/s_psi^2 .+ mkernel)

  return randn(m.K) .* sqrt.(vnew) + mnew
end
