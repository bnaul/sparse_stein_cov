library(MASS)
library(testthat)
source('sparse_stein.R')

set.seed(1)
p <- 5
n <- 1000
Sigma <- matrix(c(1.0, 0.7, 0.7, 0.0, 0.0,
                  0.7, 1.0, 0.7, 0.0, 0.0,
                  0.7, 0.7, 1.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 1.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 1.0), p, p)
X <- mvrnorm(n=n, mu=numeric(p), Sigma)
S <- t(X) %*% X

est <- sparse_stein_cov(S, n, 0.1, 1e-3)
Sigma_hat <- est$Theta
expect_equal(Sigma_hat[4:5, 1:3], matrix(0.0, 2, 3))
expect_equal(diag(Sigma_hat), diag(Sigma), tol=1e-1)
expect_equal(Sigma_hat[2:3, 1], Sigma[2:3, 1], tol=1e-1)

Omega <- solve(Sigma)
est <- sparse_stein_inv(S, n, 0.13, 1e-3)
Omega_hat <- est$Theta
expect_equal(Omega_hat[4:5, 1:3], matrix(0.0, 2, 3))
expect_equal(diag(Omega_hat), diag(Omega), tol=1e-1)
expect_equal(Omega_hat[2:3, 1], Omega[2:3, 1], tol=2e-1)
