# ----------------------------------------------------------------------
# quan_smsp: Quantile Smoothing Spline Estimator
# ----------------------------------------------------------------------
# This function estimates the alpha-th quantile function of discretely
# sampled functional data using B-splines and a roughness penalty.
# Computation is performed via fast C++ routines (IRLS + GCV).
# ----------------------------------------------------------------------

require(fda)           # For B-spline basis and penalty functions
require(Rcpp)          # Interface to C++
require(RcppArmadillo) # Efficient matrix operations in C++
Rcpp::sourceCpp("combined.cpp")  # Load the C++ functions

quan_smsp <- function(Y, alpha = 0.5, r = 2,
                      lambda_grid = exp(seq(log(1e-8), log(1e-1), length.out = 50)),
                      max_it = 100, tol = 1e-6, tun = 1e-3) {
  
  # --- Preprocessing ---
  # Convert input to matrix if needed and remove empty rows
  Y <- as.matrix(Y)
  Y <- Y[rowSums(!is.na(Y)) > 0, , drop = FALSE]
  
  n <- nrow(Y)  # number of subjects
  p <- ncol(Y)  # number of measurement points per subject
  
  # Create a uniform grid for the measurements (0 to 1)
  t_grid <- seq(0, 1, length.out = p)
  
  # Map each observation to its time point
  T_mat <- matrix(rep(t_grid, each = n), nrow = n)
  obs_idx <- which(!is.na(Y))  # indices of observed entries
  t_obs <- T_mat[obs_idx]      # time points of observed entries
  y_obs <- Y[obs_idx]          # observed values
  
  # --- Weights ---
  # Compute number of measurements per subject
  m_i <- rowSums(!is.na(Y))
  
  # Map linear indices back to row IDs
  row_id <- ((obs_idx - 1) %% n) + 1
  
  # Assign weight 1/(n * m_i) to each observation
  # Ensures subjects with fewer measurements are not overweighted
  weights_per_obs <- 1 / (n * m_i[row_id])
  
  # --- Basis Construction ---
  # Use B-spline basis with knots at observed points and order 2*r
  knots <- sort(unique(t_obs))
  b_basis <- create.bspline.basis(rangeval = c(0, 1), breaks = knots, norder = 2 * r)
  
  # Evaluate B-spline basis at observed times
  B <- eval.basis(t_obs, b_basis)
  
  # --- Penalty Matrix ---
  # Roughness penalty on r-th derivative
  Pen <- bsplinepen(b_basis, Lfdobj = r)
  
  # --- IRLS + GCV ---
  # Main C++ routine: iteratively reweighted least squares with GCV
  fit <- irls_gcv_cpp_pensp(B, Pen, y_obs, weights_per_obs, alpha,
                            lambda_grid, max_it, tol, tun)
  
  # Compute estimated quantile function on full grid
  mu_est <- eval.basis(t_grid, b_basis) %*% fit$beta_hat
  
  # --- Return Results ---
  return(list(
    mu = mu_est,       # estimated alpha-th quantile function
    lambda = fit$lambda,  # selected smoothing parameter
    weights = fit$weights # final IRLS weights
  ))
}
