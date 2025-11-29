# ----------------------------------------------------------------------
# Simulation Study: Section 4 of Kalogridis (2025+)
# ----------------------------------------------------------------------
# This script simulates discretely sampled functional data, adds
# noise, and compares three estimators:
#   1. Smoothing-spline quantile estimator (quan_smsp)
#   2. Least-squares penalized spline estimator (ls_pensp)
#   3. Quantile penalized spline estimator (quan_pensp)
# ----------------------------------------------------------------------
setwd("C:/Users/ik77w/OneDrive - University of Glasgow/Documents/GitHub/Robust-Functional-Location-Estimation")
# Number of simulations (500 in the paper)
nsim <- 50

# Number of subjects and number of measurement points
n <- 100
p <- 50

# Initialize storage for mean squared errors
mse.smsp <- mse.pensp <- mse.lspensp <- mse.huber <- rep(0, nsim)

# Initialize storage for estimated curves
shapes.smsp <- shapes.pensp <- shapes.lspensp <- shapes.huber <- matrix(0, ncol = nsim, nrow = p)

# Grid of measurement points
t_grid <- seq(0, 1, length.out = p)

# Initialize data matrix
Y <- matrix(NA, nrow = n, ncol = p)

# Define true population mean function
mu_true <- function(t) sin(2 * pi * t) # Easy example
# Harder example (mu_2 in the paper)
# mu_true <- function(t) exp(-(t-0.25)^2/0.01)+exp(-(t-0.50)^2/0.01)+exp(-(t-0.75)^2/0.01)

mu_grid <- mu_true(t_grid)

# Main Simulation Loop

for(k in 1:nsim){
  print(paste("Simulation", k, "of", nsim))
  
  # Generate functional data for n subjects
  X <- matrix(NA, nrow = n, ncol = p)
  
  for(i in 1:n){
    X[i,] <- mu_grid
    for(j in 1:50){ 
      X[i,] <- X[i,] + sqrt(2) * rt(1, df = 5) * 
        sapply(t_grid, FUN = function(x) sin((j-0.5)*pi*x)/((j-0.5)*pi))
    }
    
    # Randomly sample a subset of measurement points for subject i
    m_i <- sample(floor(0.5 * p):floor(0.8 * p), 1)
    idx <- sort(sample(seq_len(p), m_i))
    
    # Add heavy-tailed noise
    zeta <- 0.5 * rt(m_i, df = 1)
    Y[i, idx] <- X[i, idx] + zeta
  }
  # # par(mar = c(4,3.5,2,2), mgp = c(1.5, 0.75, 0), mfrow = c(1,1))
  # matplot(t_grid,t(X), pch = 20, col = "gray", cex = 1.1, xlab = "t", 
          # ylab = "", ylim = c(-4.5, 4.5), cex.lab = 1.5, cex.axis = 1.5) ; grid()
  # lines(t_grid, mu_grid, lwd = 3, col = "black")
  
  # Fit estimators
  fit.smsp <- quan_smsp(Y, alpha = 0.5)
  
  fit.lspensp <- ls_pensp(Y, K = 30)            # Least-squares penalized spline
  fit.pensp  <- quan_pensp(Y, alpha = 0.5, K = 30)  # Quantile penalized spline
  fit.huber <- huber_pensp(Y, K = 30)
  
  # -------------------------------
  # Plot a single simulation
  # -------------------------------
  # par(mar = c(4,3.5,2,2), mgp = c(3, 1.5, 0), mfrow = c(1,1))
  # plot(t_grid, mu_grid, lwd = 3, lty = 1, type = "l", cex.axis = 2.5, cex.lab = 2.5,
       # ylab = "", xlab = "t", ylim = c(-1.2, 1.2)); grid()
  # lines(t_grid, fit.huber$mu, lwd = 3, col = "blue")
  # lines(t_grid, fit.pensp$mu, lwd = 3, col = "blue")
  # lines(t_grid, fit.lspensp$mu, lwd = 3, col = "red")
  
  
  # -------------------------------
  # Compute Mean Squared Errors
  # -------------------------------
  mse.smsp[k] <- mean((fit.smsp$mu - mu_grid)^2)
  mse.pensp[k] <- mean((fit.pensp$mu - mu_grid)^2)
  mse.huber[k] <- mean((fit.huber$mu - mu_grid)^2)
  mse.lspensp[k] <- mean((fit.lspensp$mu - mu_grid)^2)
  
  # Store estimated curves
  shapes.smsp[, k] <- fit.smsp$mu
  shapes.pensp[, k] <- fit.pensp$mu
  shapes.lspensp[, k] <- fit.lspensp$mu
  shapes.huber[, k] <- fit.huber$mu
}

# Summarize Results
# Mean and standard error of MSE
mean(mse.smsp, na.rm = TRUE) * 1000; 1000 * sd(mse.smsp, na.rm = TRUE)/sqrt(nsim)
mean(mse.pensp, na.rm = TRUE) * 1000; 1000 * sd(mse.pensp, na.rm = TRUE)/sqrt(nsim)
mean(mse.lspensp, na.rm = TRUE) * 1000; 1000 * sd(mse.lspensp, na.rm = TRUE)/sqrt(nsim)
mean(mse.huber, na.rm = TRUE) * 1000; 1000 * sd(mse.huber, na.rm = TRUE)/sqrt(nsim)


# Visualize Estimated Curves Across Simulations
par(mfrow = c(1, 2))
matplot(t_grid, shapes.pensp, lwd = 3, col = "gray", type = "l", lty = 1)
lines(t_grid, mu_grid, lwd = 3, col = "black")
grid()

matplot(t_grid, shapes.lspensp, lwd = 3, col = "gray", type = "l", lty = 1)
lines(t_grid, mu_grid, lwd = 3, col = "black")
grid()

par(mfrow = c(1, 1))
matplot(t_grid, shapes.huber, lwd = 3, col = "gray", type = "l", lty = 1)
lines(t_grid, mu_grid, lwd = 3, col = "black")
grid()