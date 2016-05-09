source('util.R')

eps_cone_proj <- function(Z, eps) {
  decomp <- eigen(Z, symmetric=TRUE)
  return(decomp$vec %*% (pmax(decomp$val, eps) * t(decomp$vec)))
}

xue_cov <- function(S, lambda, eps, mu=2, tol=1e-3, max_it=500) {
  p <- nrow(S)
  Sigma <- soft_threshold(S, lambda, modify_diag=FALSE)
  Omega <- S
  Lambda <- matrix(0, p, p)
  i <- 0
  if (min(eigen(Sigma, symmetric=TRUE)$val) >= eps) {
    return(Sigma)
  }
  while (norm(Sigma - Omega) > tol) {
    if (i > max_it) {
      warning(sprintf('Maximum iterations exceeded: lambda=%f, tau=%f', lambda, eps))
      break
    }
    Omega <- eps_cone_proj(Sigma + mu * Lambda, eps)
    Sigma <- 1 / (1 + mu) * soft_threshold(mu * (S - Lambda) + Omega, lambda * mu,
                                           modify_diag=FALSE)
    Lambda <- Lambda - 1 / mu * (Omega - Sigma)
    i <- i + 1
  }
  return(Sigma)
}

liu_cov <- function(S, lambda, tau, rho=2, tol=1e-3, max_it=500) {
  p <- nrow(S)
  R <- cov2cor(S)
  Sigma <- soft_threshold(R, lambda, modify_diag=FALSE)
  Gamma <- R
  U <- matrix(0, p, p)
  i <- 0
  if (min(eigen(Sigma, symmetric=TRUE)$val) >= tau) {
    return(t(sqrt(diag(S)) * Sigma) * sqrt(diag(S)))
  }
  while (norm(Sigma - Gamma) > tol) {
    if (i > max_it) {
      warning(sprintf('Maximum iterations exceeded: lambda=%f, tau=%f', lambda, tau))
      break
    }
    Sigma <- soft_threshold(Gamma + 1 / rho * U, lambda / rho)
    Gamma <- eps_cone_proj((R + rho * Sigma - U) / (1 + rho), tau)
    U <- U + rho * (Gamma - Sigma)
    i <- i + 1
  }
  return(t(sqrt(diag(S)) * Sigma) * sqrt(diag(S)))
}
