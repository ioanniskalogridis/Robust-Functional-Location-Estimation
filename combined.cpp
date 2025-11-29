#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List ls_pensp_cpp2(const arma::mat& B,
                         const arma::mat& Pen,
                         const arma::vec& y,
                         const arma::vec& weights_per_obs,
                         Rcpp::NumericVector lambda_grid) {
  
  int m = y.n_elem;
  int k = B.n_cols;
  
  arma::vec best_beta(k, arma::fill::zeros);
  double best_gcv = arma::datum::inf;
  double best_lambda = 0;
  
  // Avoid building diagmat(W) — scale rows of B^T by weights_per_obs
  // Way more efficient...
   arma::mat Bt = B.t();
   arma::mat BtW = Bt.each_row() % weights_per_obs.t();
   arma::mat BtWB_base = BtW * B;
   arma::vec BtWy = BtW * y;

  
  for (int lg = 0; lg < lambda_grid.size(); ++lg) {
    double lambda = lambda_grid[lg];
    
    arma::mat A = BtWB_base + lambda * Pen;
    arma::vec beta = arma::solve(A, BtWy, solve_opts::fast);
    
    arma::vec fitted = B * beta;
    arma::vec resid = y - fitted;
    
    // trace_S without forming H = B * Ainv * BtW
    double trace_S = 0.0;
    for (int j = 0; j < k; ++j) {
      arma::vec col = BtWB_base.col(j);
      arma::vec x = arma::solve(A, col, solve_opts::fast);  // solve A x = column
      trace_S += x[j];                    // accumulate diagonal contribution
    }

    double denom = 1.0 - trace_S / m;
    
    double gcv = arma::mean(weights_per_obs % arma::square(resid)) / (denom * denom);
    
    if (gcv < best_gcv) {
      best_gcv = gcv;
      best_beta = beta;
      best_lambda = lambda;
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("beta") = best_beta,
    Rcpp::Named("lambda") = best_lambda,
    Rcpp::Named("gcv") = best_gcv
  );
}

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List irls_gcv_cpp_pensp(const mat& B,
                        const mat& Pen,
                        const vec& y,
                        const vec& weights_per_obs,
                        double alpha,
                        NumericVector lambda_grid,
                        int max_it,
                        double tol,
                        double tuning) {

  int m = y.n_elem;
  int k = B.n_cols;

  vec best_beta(k, fill::zeros), best_weights(m);
  double best_gcv = datum::inf, best_lambda = 0;

  // (2) Precompute B^T once
  mat B_transpose = B.t();

  for (int lg = 0; lg < lambda_grid.size(); ++lg) {
    double lambda = lambda_grid[lg];

    // Initial with only weights_per_obs
    vec beta = solve(
        (B_transpose.each_row() % weights_per_obs.t()) * B + 2 * lambda * Pen,
        B_transpose * (weights_per_obs % y),
        solve_opts::fast
    );

    vec resid = y - B * beta;
    vec weights(m);

    // Save XtW from the last IRLS iteration for the GCV
    vec w;
    mat XtW;

    int it = 0;
    for (; it < max_it; ++it) {
      // IRLS weights
      for (int i = 0; i < m; ++i) {
        if (resid[i] >= 0 && resid[i] <= tuning) {
          weights[i] = 2 * alpha / tuning;
        } else if (resid[i] < 0 && resid[i] >= -tuning) {
          weights[i] = 2 * (1 - alpha) / tuning;
        } else {
          weights[i] = (resid[i] > 0 ? alpha : alpha - 1) / resid[i];
        }
      }

      // Combined observation weights
      w = weights_per_obs % weights;

      // (1) Build XtW ONCE per iteration and reuse
      XtW = B_transpose.each_row() % w.t();

      // LHS: XtW * B
      mat XtWB = XtW * B;

      // (3) RHS without extra temp: B^T * (w % y)
      vec rhs = B_transpose * (w % y);

      vec beta_new = solve(XtWB + 2 * lambda * Pen, rhs, solve_opts::fast); // (4)
      vec resid_new = y - B * beta_new;

      if (max(abs(resid_new - resid)) < tol) {
        beta = beta_new;
        resid = resid_new;
        break;
      }

      beta = beta_new;
      resid = resid_new;
    }

    // Final GCV — reuse XtW computed in the last IRLS iteration (1)
    if (XtW.n_elem == 0) {
      vec w0 = weights_per_obs;
      XtW = B_transpose.each_row() % w0.t();
      w = w0;
    }

    mat XtWX = XtW * B + 2 * lambda * Pen;
    mat S = solve(XtWX, XtW, solve_opts::fast); // (4)
    double trace_S = trace(B * S);
    double gcv = mean(w % square(resid)) / std::pow(1 - trace_S / m, 2);

    if (gcv < best_gcv) {
      best_gcv = gcv;
      best_beta = beta;
      best_weights = weights;
      best_lambda = lambda;
    }
  }

  return List::create(
    _["beta_hat"] = best_beta,
    _["weights"] = best_weights,
    _["lambda"] = best_lambda
  );
}

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List irls_gcv_cpp_huber(const mat& B,
                        const mat& Pen,
                        const vec& y,
                        const vec& weights_per_obs,
                        NumericVector lambda_grid,
                        int max_it,
                        double tol,
                        double tuning) { 

  int m = y.n_elem;
  int k = B.n_cols;

  vec best_beta(k, fill::zeros), best_weights(m);
  double best_gcv = datum::inf, best_lambda = 0;

  mat B_transpose = B.t();

  for (int lg = 0; lg < lambda_grid.size(); ++lg) {
    double lambda = lambda_grid[lg];

    vec beta = solve(
        (B_transpose.each_row() % weights_per_obs.t()) * B + 2 * lambda * Pen,
        B_transpose * (weights_per_obs % y),
        solve_opts::fast
    );

    vec resid = y - B * beta;
    vec weights(m);

    vec w;
    mat XtW;

    int it = 0;
    for (; it < max_it; ++it) {
      // IRLS weights for Huber
      for (int i = 0; i < m; ++i) {
        if (std::abs(resid[i]) <= tuning) {
          weights[i] = 1.0;
        } else {
          weights[i] = tuning / std::abs(resid[i]);
        }
      }

      w = weights_per_obs % weights;

      XtW = B_transpose.each_row() % w.t();
      mat XtWB = XtW * B;
      vec rhs = B_transpose * (w % y);

      vec beta_new = solve(XtWB + 2 * lambda * Pen, rhs, solve_opts::fast);
      vec resid_new = y - B * beta_new;

      if (max(abs(resid_new - resid)) < tol) {
        beta = beta_new;
        resid = resid_new;
        break;
      }

      beta = beta_new;
      resid = resid_new;
    }

    if (XtW.n_elem == 0) {
      vec w0 = weights_per_obs;
      XtW = B_transpose.each_row() % w0.t();
      w = w0;
    }

    mat XtWX = XtW * B + 2 * lambda * Pen;
    mat S = solve(XtWX, XtW, solve_opts::fast);
    double trace_S = trace(B * S);
    double gcv = mean(w % square(resid)) / std::pow(1 - trace_S / m, 2);

    if (gcv < best_gcv) {
      best_gcv = gcv;
      best_beta = beta;
      best_weights = weights;
      best_lambda = lambda;
    }
  }

  return List::create(
    _["beta_hat"] = best_beta,
    _["weights"] = best_weights,
    _["lambda"] = best_lambda
  );
}
