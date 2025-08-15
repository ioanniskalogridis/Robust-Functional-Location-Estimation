require(fda) 
require(Rcpp) 
require(RcppArmadillo)
Rcpp::sourceCpp("combined.cpp")
ls_pensp <- function(Y, r = 2, m = 4, K = 30, lambda_grid = exp(seq(log(1e-8), log(1e-1), 
                                                                    length.out = 50))) { 
  
  Y <- as.matrix(Y) 
  Y <- Y[rowSums(!is.na(Y)) > 0, , drop = FALSE] 
  n <- nrow(Y) 
  p <- ncol(Y)
  
  t_grid <- seq(0, 1, length.out = p) 
  T_mat <- matrix(rep(t_grid, each = n), nrow = n) 
  obs_idx <- which(!is.na(Y)) 
  t_obs <- T_mat[obs_idx] 
  y_obs <- as.numeric(Y[obs_idx])
  
  m_i <- rowSums(!is.na(Y)) 
  row_id <- ((obs_idx - 1) %% n) + 1 
  weights_per_obs <- 1 / (n * m_i[row_id])
  
  b_basis <- create.bspline.basis(c(0, 1), nbasis = K, norder = m) 
  B <- eval.basis(t_obs, b_basis) 
  Pen <- bsplinepen(b_basis, Lfdobj = r)
  
  res <- ls_pensp_cpp2(B, Pen, y_obs, weights_per_obs, lambda_grid)
  
  mu_est <- as.numeric(eval.basis(t_grid, b_basis) %*% res$beta)
  
  list(mu = mu_est, beta = res$beta, lambda = res$lambda, 
       gcv = res$gcv, t_grid = t_grid, basis = b_basis, Penalty = Pen, 
       weights_per_obs = weights_per_obs) }