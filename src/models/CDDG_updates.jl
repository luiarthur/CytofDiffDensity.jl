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

function update_lambdaC(m::CDDG, s::T) where T
  # TODO
  return s.lambdaC
end

function update_lambdaT(m::CDDG, s::T) where T
  # TODO
  return s.lambdaT
end

function update_vC(m::CDDG, s::T) where T
  # TODO
  return s.vC
end

function update_vT(m::CDDG, s::T) where T
  # TODO
  return s.vT
end

function update_zetaC(m::CDDG, s::T) where T
  # TODO
  return s.zetaC
end

function update_zetaT(m::CDDG, s::T) where T
  # TODO
  return s.zetaT
end

function update_mu(m::CDDG, s::T) where T
  # TODO
  return s.mu
end

function update_tau(m::CDDG, s::T) where T
  # TODO
  return s.tau
end

function update_omega(m::CDDG, s::T) where T
  # TODO
  return s.omega
end

function logprob_nu(m::CDDG, s::T, nu::AbstractVector{<:Real}) where T
  # TODO
  logprior = sum(logpdf.(m.nu, nu))
  loglike = 0
  return logprior + loglike
end

function update_psi(m::CDDG, s::T) where T
  # TODO
  return s.psi
end
