# ----------------------------------------------------------------------
# Simulation Study: Section 4 of Kalogridis (2025+)
# ----------------------------------------------------------------------
# This script simulates discretely sampled functional data, adds
# noise, and compares 4 estimators:
#   1. LAD smoothing spline estimator (quan_smsp)
#   2. Least-squares penalized spline estimator (ls_pensp)
#   3. LAD penalized spline estimator (quan_pensp)
#   4. Huber penalized spline estimator (huber_pensp)
# ----------------------------------------------------------------------
# Number of simulations (500 in the paper)
nsim <- 5

# Number of subjects and number of measurement points
n <- 100 # 50
p <- 50 # 30

# Initialize storage for mean squared errors
mse.smsp <- mse.pensp <- mse.lspensp <- mse.huber <- rep(0, nsim)

# Initialize storage for estimated curves
shapes.smsp <- shapes.pensp <- shapes.lspensp <- shapes.huber <- matrix(0, ncol = nsim, nrow = p)
shapes.sfpca <- matrix(0, ncol = nsim, nrow = p)

# Grid of measurement points
t_grid <- seq(0, 1, length.out = p)

# True population mean functions
mu_true <- function(t) sin(2 * pi * t) # Easy example
# Harder example (mu_2 in the paper)
# mu_true <- function(t) exp(-(t-0.25)^2/0.01)+exp(-(t-0.50)^2/0.01)+exp(-(t-0.75)^2/0.01)
mu_grid <- mu_true(t_grid)# Evaluate at grid

# Main Simulation Loop

for(k in 1:nsim){
  print(paste("Simulation", k, "of", nsim))
  
  # Generate functional data for n subjects
  X <- matrix(NA, nrow = n, ncol = p) # discretely sampled data
  Y <- matrix(NA, nrow = n, ncol = p) # discretely sampled data with measurement error
  
  for(i in 1:n){
    X[i,] <- mu_grid
    for(j in 1:50){ 
      X[i,] <- X[i,] + sqrt(2) * rt(1, df = 5) * # or df = 1 for very heavy-tailed curves
        sapply(t_grid, FUN = function(x) sin((j-0.5)*pi*x)/((j-0.5)*pi))
    }
    
    # Randomly sample a subset of measurement points for subject i
    m_i <- sample(floor(0.5 * p):floor(0.8 * p), 1)
    idx <- sort(sample(seq_len(p), m_i))
    
    # Add heavy-tailed noise or set zeta = 0
    zeta <- 0.5 * rt(m_i, df = 5)
    # zeta <- 0
    Y[i, idx] <- X[i, idx] + zeta
  }
  # Fit estimators
  fit.smsp <- quan_smsp(Y) # LAD smoothing splines
  fit.lspensp <- ls_pensp(Y)            # Least squares penalized splines
  fit.pensp  <- quan_pensp(Y)  # LAD penalized splines
  fit.huber <- huber_pensp(Y) # Huber penalized splines
  
  # Plot a single simulation
  par(mar = c(4,3.5,2,2), mgp = c(3, 1.5, 0), mfrow = c(1,1))
  plot(t_grid, mu_grid, lwd = 3, lty = 1, type = "l", cex.axis = 2.5, cex.lab = 2.5,
       ylab = "", xlab = "t", ylim = c(-1.2, 1.2)); grid()
  lines(t_grid, fit.huber$mu, lwd = 3, col = "blue")
  lines(t_grid, fit.pensp$mu, lwd = 3, col = "orange")
  lines(t_grid, fit.lspensp$mu, lwd = 3, col = "red")
  
  
  # Compute Mean Squared Errors
  mse.smsp[k] <- mean((fit.smsp$mu - mu_grid)^2)
  mse.pensp[k] <- mean((fit.pensp$mu - mu_grid)^2)
  mse.lspensp[k] <- mean((fit.lspensp$mu - mu_grid)^2)
  mse.huber[k] <- mean((fit.huber$mu - mu_grid)^2)
  
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
# 
par(mar = c(4,3.5,2,2), mgp = c(1.5, 0.75, 0), mfrow = c(1,1))
matplot(t_grid, shapes.pensp, lwd = 3, col = "gray", type = "l", lty = 1, cex.lab = 1.5, cex.axis = 1.5,
        ylab = "", xlab = "t", ylim = c(-1,1)) # ylim = c(-1,1) or ylim = c(-0.5,1.3)
lines(t_grid, mu_grid, lwd = 3, col = "black")
grid()

par(mar = c(4,3.5,2,2), mgp = c(1.5, 0.75, 0), mfrow = c(1,1))
matplot(t_grid, shapes.lspensp, lwd = 3, col = "gray", type = "l", lty = 1,cex.lab = 1.5, cex.axis = 1.5,
        ylab = "", xlab = "t", ylim = c(-1, 1))
lines(t_grid, mu_grid, lwd = 3, col = "black")
grid()

par(mfrow = c(1, 1))
matplot(t_grid, shapes.huber, lwd = 3, col = "gray", type = "l", lty = 1)
lines(t_grid, mu_grid, lwd = 3, col = "black")
grid()