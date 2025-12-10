# ----------------------------------------------------------------------
# quan_pensp: Quantile Penalized Spline Estimator
# ----------------------------------------------------------------------
# This function estimates the conditional quantile function of 
# discretely sampled functional data using B-splines and an O-spline roughness 
# penalty. Computation is performed via a fast C++ routine.
# For details, see the documentation below.
# ----------------------------------------------------------------------

require(fda)           # For B-spline basis and penalty functions
require(Rcpp)          # Interface to C++
require(RcppArmadillo) # Efficient matrix operations in C++
Rcpp::sourceCpp("combined.cpp")  # Load the C++ functions

quan_pensp <- function(Y, alpha = 0.5, r = 2, m = 4, K = NULL,
                       lambda_grid = exp(seq(log(1e-8), log(1e-1), length.out = 50)),
                       max_it = 200, tol = 1e-6, tun = 1e-3) {
  
  # Y: Numeric matrix (subjects x time points). NA for missing values.
  # alpha: quantile to be estimated (default = 0.5, median)
  # r: Order of the penalty (default = 2, 2nd derivative penalization)
  # m: Order of B-spline basis (default = 4, cubic splines)
  # K: Number of interior knots (default = min(35, max observed per subject))
  # lambda_grid: Candidate penalty parameters for GCV selection
  # max_it: Maximum IRLS iterations (default = 200)
  # tol: Numeric tolerance for IRLS convergence (default = 1e-6)
  # tun: tuning of the local quadratic approximation to the check loss (default = 10^{-3})
  
  # - Preprocessing -
  # Convert input to matrix and remove rows with all NA
  Y <- as.matrix(Y)
  Y <- Y[rowSums(!is.na(Y)) > 0, , drop = FALSE]
  
  n <- nrow(Y)  # number of subjects
  p <- ncol(Y)  # maximum number of measurement points 
  
  # Uniform grid of measurement points from 0 to 1
  t_grid <- seq(0, 1, length.out = p)
  
  # Map observations to their time points
  T_mat <- matrix(rep(t_grid, each = n), nrow = n)
  obs_idx <- which(!is.na(Y))      # indices of observed entries
  t_obs <- T_mat[obs_idx]          # time points of observed entries
  y_obs <- Y[obs_idx]              # observed values

  # Number of measurements per subject
  m_i <- rowSums(!is.na(Y))
  K <- ifelse(is.null(K), min(35, max(m_i)), K)
          
  # Map linear indices back to row IDs
  row_id <- ((obs_idx - 1) %% n) + 1
  
  # Assign weight 1/(n * m_i) to each observation
  weights_per_obs <- 1 / (n * m_i[row_id])
  
  # Create B-spline basis with K basis functions and order m
  b_basis <- create.bspline.basis(rangeval = c(0, 1), nbasis = K, norder = m)
  
  # Evaluate B-spline basis at observed points
  B <- eval.basis(t_obs, b_basis)
  
  # Penalty Matrix from the fda package 
  # Roughness penalty on r-th derivative
  Pen <- bsplinepen(b_basis, Lfdobj = r)
  
  # IRLS + GCV
  # Fit the penalized quantile regression using the C++ routine
  fit <- irls_gcv_cpp_pensp(B, Pen, y_obs, weights_per_obs, alpha,
                            lambda_grid, max_it, tol, tun)
  
  # Estimated quantile function on full grid
  mu_est <- eval.basis(t_grid, b_basis) %*% fit$beta_hat
  
  return(list(
    mu = mu_est,       # estimated quantile function
    lambda = fit$lambda,  # selected smoothing parameter
    weights = fit$weights, # converged IRLS weights
    iter = fit$iterations # number of iterations for the best lambda
  ))
}
