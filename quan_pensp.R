require(fda) 
require(Rcpp) 
require(RcppArmadillo)
setwd("C:/Users/ik77w/Documents")
Rcpp::sourceCpp("irls_gcv_cpp_pensp.cpp")
# Rcpp::sourceCpp("combined.cpp")
quan_pensp <- function(Y, alpha = 0.5, r = 2, m = 4, K = 30,
                       lambda_grid = exp(seq(log(1e-8), log(1e-1), length.out = 50)),
                       max_it = 100, tol = 1e-6, tun = 1e-3) {
  
  # --- Preprocessing ---
  Y <- as.matrix(Y)
  Y <- Y[rowSums(!is.na(Y)) > 0, , drop = FALSE]
  n <- nrow(Y)
  p <- ncol(Y)
  
  t_grid <- seq(0, 1, length.out = p)
  T_mat <- matrix(rep(t_grid, each = n), nrow = n)
  obs_idx <- which(!is.na(Y))
  t_obs <- T_mat[obs_idx]
  y_obs <- Y[obs_idx]
  
  # --- Weights ---
  m_i <- rowSums(!is.na(Y))  # number of obs per subject
  row_id <- ((obs_idx - 1) %% n) + 1  # maps linear index to row
  weights_per_obs <- 1 / (n * m_i[row_id])
  
  # --- Basis Construction ---
  knots <- sort(unique(t_obs))
  b_basis <- create.bspline.basis(rangeval = c(0, 1), nbasis = K, norder = m)
  B <- eval.basis(t_obs, b_basis)
  
  # --- Penalty ---
  Pen <- bsplinepen(b_basis, Lfdobj = r)
  
  # --- Initial Fit ---
  h_m <- 1 / mean(1 / m_i)
  lambda_init <- 100 * (h_m * n)^(-2 * r / (2 * r + 1))
  beta_init <- solve(t(B) %*% diag(weights_per_obs) %*% B + lambda_init * Pen,
                     t(B) %*% (weights_per_obs * y_obs))
  resid_init <- y_obs - B %*% beta_init
  
  # --- IRLS + GCV ---
  fit <- irls_gcv_cpp_pensp(B, Pen, y_obs, weights_per_obs, alpha,
                            lambda_grid, max_it, tol, tun)
  
  mu_est <- eval.basis(t_grid, b_basis) %*% fit$beta_hat
  
  return(list(mu = mu_est,
              lambda = fit$lambda,
              weights = fit$weights))
}
