# Sparse Steinian Covariance Estimation
- `sparse_stein.R` implements the sparse covariance estimator of Naul and Taylor, *Sparse
  Steinian Covariance Estimation*, Journal of Computational Graphical Statistics (2016).
    - `sparse_stein_cov(...)` returns the sparse Stein estimate for the covariance matrix
    - `sparse_stein_inv(...)` returns the sparse Stein estimate for the concentration matrix
    - Example: 

       ```R
       p <- 6
       n <- 15
       S <- rWishart(1, n, diag(p))[,,1]
       Sigma_hat <- sparse_stein_cov(S, n, lambda=0.5, rho=0.5)$Theta
       ```
- `sparse_mle.R` implements two other sparse covariance estimators for comparison, namely: 
    - Xue et al., *Positive-definite l1-penalized estimation of large covariance matrices*, Journal of the American Statistical Association (2012).
    - Liu et al., *Sparse covariance matrix estimation with eigenvalue constraints*, Journal of Computational and Graphical Statistics (2014).
- `simulations.R` contains code used to produce the figures in Naul and Taylor, *Sparse
  Steinian Covariance Estimation*, Journal of Computational Graphical Statistics (2016). See comments in that file for details.
