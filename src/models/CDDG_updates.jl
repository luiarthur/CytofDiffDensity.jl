function grab(from::T, p::Symbol, sample_idx::Symbol) where T
  return getfield(from, Symbol(p, sample_idx))
end

function update_etaC(m::CDDG, s::T) where T
  anew = copy(m.eta.alpha)

  foreach(lam -> anew[lam] +=1, s.lambdaC)
  m.beta || foreach(lam -> anew[lam] +=1, s.lambdaT)

  return rand(Dirichlet(anew))
end

function update_etaT(m::CDDG, s::T) where T
  # FIXME: The manuscript is wrong. Need to correct.
  if m.beta
    anew = copy(m.eta.alpha)
    foreach(lam -> anew[lam] +=1, s.lambdaT)
    return rand(Dirichlet(anew))
  else
    return m.etaC
  end
end

function update_lambda(m::CDDG, s::T, sample_idx::Symbol) where T
  y = grab(m, :y, sample_idx)
  v = grab(s, :v, sample_idx)
  zeta = grab(s, :zeta, sample_idx)
  eta = grab(s, :eta, sample_idx)
  logeta = log.(eta)

  lambda = [let
              loc = s.mu + s.psi * zeta[n]
              scale = sqrt.(s.omega) / sqrt(v[n])
              logmix = normlogpdf.(loc, scale, y[n]) + logeta
              all(s.nu .< Inf) && (logmix += gammalogpdf.(s.nu/2, 2 ./ s.nu, v[n]))
              MCMC.wsample_logprob(logmix)
            end for n in eachindex(y)]

  return lambda
end

update_lambdaC(m::CDDG, s::T) where T = update_lambda(m, s, :C)
update_lambdaT(m::CDDG, s::T) where T = update_lambda(m, s, :T)

function update_v(m::CDDG, s::T, i::Symbol) where T
  lami = grab(s, :lambda, i)
  zetai = grab(s, :zeta, i)
  yi = grab(m, :y, i)
  vi = grab(s, :v, i)

  nu = s.nu[lami]
  omega = s.omega[lami]
  psi = s.psi[lami]
  mu = s.mu[lami]

  shape = nu/2 .+ 1
  rate = (nu + zetai.^2 + ((yi - mu - psi .* zetai) .^ 2) ./ omega) / 2

  return @.rand(Gamma(shape, 1 / rate))
end

update_vC(m::CDDG, s::T) where T = update_v(m, s, :C)
update_vT(m::CDDG, s::T) where T = update_v(m, s, :T)

function update_zeta(m::CDDG, s::T, i::Symbol) where T
  yi = grab(m, :y, i)
  zetai = grab(s, :zeta, i)
  lami = grab(s, :lambda, i)
  vi = grab(s, :v, i)

  psi = s.psi[lami]
  omega = s.omega[lami]
  mu = s.mu[lami]

  vnew = @. 1 / (vi + (psi ^ 2) * vi / omega)
  mnew = @. vnew * vi * psi * (yi - mu) / omega

  return @.rand(truncated(Normal(mnew, sqrt(vnew)), 0, Inf))
end

update_zetaC(m::CDDG, s::T) where T = update_zeta(m, s, :C)
update_zetaT(m::CDDG, s::T) where T = update_zeta(m, s, :T)

function update_mu(m::CDDG, s::T) where T
  m_mu, s_mu = params(m.mu)
  vkernels = zero.(s.mu)
  mkernels = zero.(s.mu)

  for i in (:C, :T)
    y = grab(m, :y, i)
    lambda = grab(s, :lambda, i)
    zeta = grab(s, :zeta, i)
    v = grab(s, :v, i)
    for n in eachindex(lambda)
      k = lambda[n]
      vkernels[k] += v[n]
      mkernels[k] += (y[n] - s.psi[k] * zeta[n]) * v[n]
    end
  end

  vnew = 1 ./ (s_mu^-2 .+ vkernels ./ s.omega)
  mnew = vnew .* (m_mu/(s_mu^2) .+ mkernels ./ s.omega)
  
  return randn(m.K) .* sqrt.(vnew) + mnew
end

function update_tau(m::CDDG, s::T) where T
  a, b = params(m.tau)
  new_shape = a + m.K * m.a_omega
  new_rate = b + sum(1 ./ s.omega)
  return rand(Gamma(new_shape, 1 / new_rate))
end

function update_sigma!(s::T) where T
  @. s.sigma = scalefromaltskewt(sqrt(s.omega), s.psi)
end

function update_phi!(s::T) where T
  @. s.phi = skewfromaltskewt(sqrt(s.omega), s.psi)
end

function update_sigma_phi!(s::T) where T
  update_sigma!(s)
  update_phi!(s)
end

function update_omega(m::CDDG, s::T) where T
  a, b = m.a_omega, s.tau
  akernel = zero.(s.omega)
  bkernel = zero.(s.omega)

  for i in (:C, :T)
    y = grab(m, :y, i)
    lam = grab(s, :lambda, i)
    zeta = grab(s, :zeta, i)
    v = grab(s, :v, i)
    for n in eachindex(y)
      k = lam[n]
      akernel[k] += 1
      bkernel[k] += v[n] * (y[n] - s.mu[k] - s.psi[k] * zeta[n])^2
    end
  end

  anew = a .+ akernel / 2
  bnew = b .+ bkernel / 2

  @. s.omega = rand(InverseGamma(anew, bnew))

  update_sigma_phi!(s)

  return s.omega
end

function logprob_nu(m::CDDG, s::T, nu::AbstractVector{<:Real}) where T
  logprior = sum(logpdf.(m.nu, nu))
  loglike = let
    _ll = 0.0
    for i in (:C, :T)
      lam = grab(s, :lambda, i)
      v = grab(s, :v, i)
      nu_i = nu[lam]
      _ll += sum(gammalogpdf.(nu_i / 2, 2 ./ nu_i, v))
    end
    _ll
  end
  return logprior + loglike
end

function update_psi(m::CDDG, s::T) where T
  m_psi, s_psi = params(m.psi)
  vkernel = zero.(s.psi)
  mkernel = zero.(s.psi)

  # Update kernels
  for i in (:C, :T)
    yfinite = grab(m, :y, i)
    lam = grab(s, :lambda, i)
    zeta = grab(s, :zeta, i)
    v = grab(s, :v, i)
    for n in eachindex(v)
      k = lam[n]
      vkernel[k] += zeta[n]^2 * v[n] / s.omega[k]
      mkernel[k] += zeta[n] * (yfinite[n] - s.mu[k]) * v[n] / s.omega[k]
    end
  end

  vnew = 1 ./ (s_psi^-2 .+ vkernel)
  mnew = vnew .* (m_psi/s_psi^2 .+ mkernel)

  s.psi .= randn(m.K) .* sqrt.(vnew) + mnew

  update_sigma_phi!(s)

  return s.psi
end
