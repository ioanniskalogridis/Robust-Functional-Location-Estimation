## --- Simulation with heavy-tailed KL process --- ##
nsim <- 500
n    <- 100
p    <- 50

mse.smsp <- mse.pensp <- mse.lspensp <- rep(0, nsim)
shapes.smsp <- shapes.pensp <- shapes.lspensp <- matrix(0, ncol = nsim, nrow = p)

t_grid <- seq(0, 1, length.out = p)
Y <- matrix(NA, nrow = n, ncol = p)

## Population mean
mu_true <- function(t) sin(2 * pi * t) # Easy
# mu_true <- function(t) exp(-(t-0.25)^2/0.01)+exp(-(t-0.50)^2/0.01) + exp(-(t-0.75)^2/0.01) # Much harder

mu_grid <- mu_true(t_grid)

for(k in 1:nsim){
  print(k)
  X <- matrix(NA, nrow = n, ncol = p)
  for(i in 1:n){
    X[i,] <- mu_grid 
    for(j in 1:50){ 
      X[i,] <- X[i, ] +  sqrt(2)*rt(1, df = 5)*sapply(t_grid, FUN = function(x) sin((j-1/2)*pi*x)/((j-1/2)*pi) )
      # X[i,] <- X[i, ] +  sqrt(2)*rnorm(1)*sapply(t_grid, FUN = function(x) sin((j-1/2)*pi*x)/((j-1/2)*pi) )
    }
    m_i <- sample(floor(0.5 * p):floor(0.8 * p), 1)
    idx <- sort(sample(seq_len(p), m_i))
    zeta <- 0.5*rt(m_i, df = 10)
    Y[i, idx] <- X[i, idx] + zeta
  }
  # par(mar = c(4,3.5,2,2), mgp = c(3, 1.5, 0))
  # matplot(t_grid,t(Y), lwd = 3, lty = 1, col = "gray", type = "p", pch = 20,
          # cex.lab = 2.5, cex.axis = 2.5, ylab = "", xlab = "t", cex = 1.6)
  # lines(t_grid, mu_grid, lwd = 3, col = "black")
  # grid()

  fit.smsp  <- quan_smsp(Y,  alpha = 0.5)         #  smoothing-spline quantile estimator
  fit.lspensp <- ls_pensp(Y, K = 30)
  fit.pensp <- quan_pensp(Y, alpha = 0.5, K = 30) # penalized-spline quantile estimator
  
  plot()
  plot(t_grid,fit.pensp$mu, lwd = 3, type= "l", col = "blue")
  lines(t_grid, fit.lspensp$mu, lwd = 3, type = "l", col = "red")
  
  ## --- Compute MSEs against population mean mu_grid ---
  # If you want to compare to the population mean mu_grid:
  mse.smsp[k]  <- mean((fit.smsp$mu - mu_grid)^2)
  mse.pensp[k] <- mean((fit.pensp$mu - mu_grid)^2)
  mse.lspensp[k] <- mean((fit.lspensp$mu - mu_grid)^2)
  
  ## --- Save fitted shapes for later plotting ---
  shapes.smsp[, k]  <- fit.smsp$mu
  shapes.pensp[, k] <- fit.pensp$mu
  shapes.lspensp[, k] <- fit.lspensp$mu
}

mean(mse.smsp, na.rm = TRUE)*1000; 1000*sd(mse.smsp, na.rm = TRUE)/sqrt(nsim)
mean(mse.pensp, na.rm = TRUE)*1000;1000*sd(mse.pensp, na.rm = TRUE)/sqrt(nsim)
mean(mse.lspensp, na.rm = TRUE)*1000;1000*sd(mse.lspensp, na.rm = TRUE)/sqrt(nsim)

matplot(t_grid,shapes.smsp,lwd = 3, col = "gray", type = "l", lty = 1); lines(t_grid, mu_grid, lwd = 3, col = "black")
matplot(t_grid,shapes.pensp,lwd = 3, col = "gray", type = "l", lty = 1); lines(t_grid, mu_grid, lwd = 3, col = "black")
matplot(t_grid,shapes.lspensp,lwd = 3, col = "gray", type = "l", lty = 1); lines(t_grid, mu_grid, lwd = 3, col = "black")

par(mar = c(4,3.5,2,2), mgp = c(3, 1.5, 0))
# matplot(t_grid,shapes.smsp,lwd = 3, col = "gray", type = "l", lty = 1); lines(t_grid, mu_grid, lwd = 3, col = "black")
matplot(t_grid,shapes.pensp,lwd = 3, col = "gray", type = "l", 
        lty = 1, cex.lab = 2.5, cex.axis = 2.5, ylab = "", xlab = "t"); lines(t_grid, mu_grid, lwd = 3, col = "black")
grid()
matplot(t_grid,shapes.lspensp,lwd = 3, col = "gray", type = "l", lty = 1,
        cex.lab = 2.5, cex.axis = 2.5, ylab = "", xlab = "t"); lines(t_grid, mu_grid, lwd = 3, col = "black")
grid()
