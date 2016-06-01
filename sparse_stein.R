library(Iso)
source('util.R')

# vector with ith component \sum_{j \neq i} 1/(l_i-l_j); appears in UBEOR expression
sum_diffs <- function(l) {rowSums(1/(outer(l, l, '-') + diag(Inf, length(l))))}
# vector with ith component \sum_{j \neq i} l_j/(l_i-l_j)^2; appears in UBEOR expression
sum_diffs_lj <- function(l) {rowSums(t(l/(t(outer(l, l, '-')^2) + diag(Inf, length(l)))))}

# Computes the sparse-Stein covariance estimate for the given sample covariance
# Parameters
# ----------
# S: sample covariance matrix (UNSCALED, i.e. X'T %*% X)
# n: number of samples
# lambda: l_1 penalty parameter
# rho: Frobenius norm penalty parameter
# delta: "adaptive lasso" penalty parameter; increase to penalize small values
#        more strongly ("adaptive Sparse Stein" <--> delta=2)
# Outputs
# -------
# Theta: estimated covariance
# phi: eigenvalues of un-thresholded estimate
# H: eigenvectors of un-thresholded estimate
sparse_stein_cov <- function(S, n, lambda, rho, delta=0, tol=1e-3, max_it=500) {
  p <- nrow(S)
  decomp <- eigen(S, symmetric=TRUE)
  l <- decomp$val; H <- decomp$vec
  # only use non-zero eigenvalues/corresponding eigenvectors
  if (p > n) {
    l <- l[1:n]; H <- H[,1:n]
  }

  # A = H Phi H^T (dense estimate); intialize @ Stein/Haff's estimator
  # Theta = soft_threshold(A); initialize @ thresholded Stein/Haff's
  alpha <- abs(n - p) + 1 + 2 * l * sum_diffs(l)
  phi <- 1 / pava(alpha / l)
  A <- H %*% (phi * t(H))
  Lambda <- lambda / median(abs(A) ^ -delta) * abs(A) ^ (-delta)
  Theta <- soft_threshold(A, Lambda)

  i <- 0
  if (rho > 0) {  # for rho=0, optimal value is just thresholded Stein
    repeat {
      if (i > max_it) {
        warning(sprintf('Maximum iterations exceeded: lambda=%f, rho=%f', lambda, rho))
        break
      }
      i <- i+1
      x <- rowSums((t(H) %*% Theta) * t(H)) # efficient version of diag(t(H) %*% Theta %*% H)
      phi_old <- phi
      ystar <- -pava(-(rho * x - alpha / l))
      phi <- 1 / 2 * (ystar / rho + sqrt(ystar ^ 2 / rho ^ 2 + 4 / rho))
      A <- H %*% (phi * t(H))
      Theta <- soft_threshold(A, Lambda)
      if (sqrt(sum((phi - phi_old) ^ 2)) / p < tol) {
        break
      }
    }
  }
  return(list(phi=phi, H=H, Theta=Theta))
}

# Computes the sparse-Stein inverse covariance estimate for the given sample covariance
# Parameters
# ----------
# S: sample covariance matrix (UNSCALED, i.e. X'T %*% X)
# n: number of samples
# lambda: l_1 penalty parameter
# rho: Frobenius norm penalty parameter
# delta: "adaptive lasso" penalty parameter; increase to penalize small values
#        more strongly ("adaptive Sparse Stein" <--> delta=2)
# Outputs
# -------
# Theta: estimated inverse covariance
# phi: eigenvalues of un-thresholded estimate
# H: eigenvectors of un-thresholded estimate
sparse_stein_inv <- function(S, n, lambda, rho, loss='F', delta=0, tol=1e-3, max_it=500) {
  p <- nrow(S)
  decomp <- eigen(S, symmetric=TRUE)
  l <- decomp$val; H <- decomp$vec
  # only use non-zero eigenvalues/corresponding eigenvectors
  if (p > n) {
    l <- l[1:n]; H <- H[,1:n]
  }

  # A = H Phi H^T (dense estimate); intialize @ Stein/Haff's estimator
  # Theta = soft_threshold(A); initialize @ thresholded Stein/Haff's
  if (loss == 'F' || loss == 'L0') {
    alpha <- abs(n - p) - 3 + 2 * l * sum_diffs(l)
    g <- alpha / l
    w <- l ^ 0
  } else if (loss == 'L1') {
    alpha <- abs(n - p) - 1 + 2 * l * sum_diffs(l)
    g <- alpha / l
    w <- l ^ 1
  } else if (loss == 'L2') {
    alpha <- abs(n - p) + 1 + 2 * l * sum_diffs(l)
    g <- alpha / l
    w <- l ^ 2
  } else stop('Unrecognized loss function')
  phi <- pava(g, w)
  A <- H %*% diag(phi) %*% t(H)
  Lambda <- lambda / median(abs(A) ^ -delta) * rho * (abs(A) ^ -delta)
  Theta <- soft_threshold(A, Lambda / rho)

  i <- 0
  if (rho > 0) {  # for rho=0, optimal value is just thresholded Stein
    repeat {
      if (i > max_it) {
        warning(sprintf('Maximum iterations exceeded: lambda=%f, rho=%f', lambda, rho))
        break
      }
      i <- i+1
      x <- rowSums((t(H) %*% Theta) * t(H)) # efficient version of diag(t(H) %*% Theta %*% H)
      phi_old <- phi
      phi <- pava(1 / (w + rho / 2) * (g * w + rho / 2 * x), w + rho / 2)
      A <- H %*% diag(phi) %*% t(H)
      Theta <- soft_threshold(A, Lambda / rho)
      if (sqrt(sum((phi - phi_old)^2)) < tol) {
        break
      }
    }
  }
  return(list(phi=phi, H=H, Theta=Theta))
}
