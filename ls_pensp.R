# ----------------------------------------------------------------------
# ls_pensp: Least-Squares Penalized Spline Estimator
# ----------------------------------------------------------------------
# This function estimates the mean function of discretely sampled
# functional data using B-splines and an O-spline roughness penalty.
# Computation is performed via a fast C++ routine.
# For details, see the documentation below.
# ----------------------------------------------------------------------

require(fda)           # For B-spline basis and penalty functions
require(Rcpp)          # Interface to C++
require(RcppArmadillo) # Efficient matrix operations in C++
Rcpp::sourceCpp("combined.cpp")  # Load the C++ functions

ls_pensp <- function(Y, r = 2, m = 4, K = NULL,
                     lambda_grid = exp(seq(log(1e-8), log(1e-1), length.out = 50))) { 
  
  # Y: Numeric matrix (subjects x time points). NA for missing values.
  # r: Order of the penalty (default = 2, 2nd derivative penalization)
  # m: Order of B-spline basis (default = 4, cubic splines)
  # K: Number of interior knots (default = min(35, max observed per subject))
  # lambda_grid: Candidate penalty parameters for GCV selection
  # max_it: Maximum IRLS iterations (default = 200)
  # tol: Numeric tolerance for IRLS convergence (default = 1e-6)
  
  # --- Preprocessing ---
  # Convert input to matrix if needed and remove empty rows
  Y <- as.matrix(Y)
  Y <- Y[rowSums(!is.na(Y)) > 0, , drop = FALSE]
  
  n <- nrow(Y)  # number of subjects
  p <- ncol(Y)  # number of measurement points per subject
  
  # Create a uniform grid for the measurements (0 to 1)
  t_grid <- seq(0, 1, length.out = p)
  
  # Map each observation to its discretization point
  T_mat <- matrix(rep(t_grid, each = n), nrow = n)
  obs_idx <- which(!is.na(Y))      # indices of observed entries
  t_obs <- T_mat[obs_idx]          # time points of observed entries
  y_obs <- as.numeric(Y[obs_idx])  # observed values

  # --- Weights
  # Compute number of measurements per subject
  m_i <- rowSums(!is.na(Y))
  K <- ifelse(is.null(K), min(35, max(m_i)), K)
  
  # Map linear indices back to row IDs
  row_id <- ((obs_idx - 1) %% n) + 1
  
  # Assign weight 1/(n * m_i) to each observation
  weights_per_obs <- 1 / (n * m_i[row_id])
  
  # Use B-spline basis with K basis functions from equidistant knots and order m
  b_basis <- create.bspline.basis(rangeval = c(0, 1), nbasis = K, norder = m)
  
  # Evaluate B-spline basis at observed times
  B <- eval.basis(t_obs, b_basis)
  
  # --- Penalty Matrix ---
  # Roughness penalty on r-th derivative
  Pen <- bsplinepen(b_basis, Lfdobj = r)
  
  # Solve penalized least squares using the fast C++ routine
  res <- ls_pensp_cpp2(B, Pen, y_obs, weights_per_obs, lambda_grid)
  
  # Compute estimated mean function on full grid
  mu_est <- as.numeric(eval.basis(t_grid, b_basis) %*% res$beta)
  
  # --- Return Results ---
  return(list(
    mu = mu_est,                  # estimated mean function
    beta = res$beta,              # estimated B-spline coefficients
    lambda = res$lambda,          # selected smoothing parameter
    gcv = res$gcv,                # GCV values over lambda grid
    t_grid = t_grid,              # grid of evaluation points
    basis = b_basis,              # B-spline basis object
    Penalty = Pen,                # penalty matrix
    weights_per_obs = weights_per_obs  # observation weights
  ))
}
