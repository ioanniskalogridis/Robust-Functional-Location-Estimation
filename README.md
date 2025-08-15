This repository contains fast C++ implementations with an R interface of the penalized spline estimators of Kalogridis (2025+) as well as a remake of the older smoothing spline estimator of
Kalogridis and Van Aelst (2023, SJS). 

Here are detailed instructions:
1. First, download all the files in your R working directory.
2. Open and load the R-functions quad_smsp.R (Quantile Smoothing Spline Estimator), quad_pensp.R (Quantile Penalized Spline Estimator) and ls_pensp2.R (Least-Squares Penalized Spline Estimator).
3. Be sure to have installed and loaded the R-packages fda, Rcpp and RcppArmadillo.
4. These R functions will source the combined.cpp file containing the C++ implementations; no C++ knowledge is required.

Here is an example of their use in R with discretely sampled functional data:

