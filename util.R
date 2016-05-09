soft_threshold <- function(A, eps, modify_diag=FALSE, scaled=FALSE) {
  if (scaled) {A_eps <- cov2cor(A)} else {A_eps <- A}
  A_eps <- sign(A_eps) * pmax(abs(A_eps) - eps, 0)
  if (scaled) A_eps <- t(A_eps * sqrt(diag(A))) * sqrt(diag(A))
  if (!modify_diag) diag(A_eps) <- diag(A)
  return(A_eps)
}
