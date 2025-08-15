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
  
  # plot(t_grid,fit.pensp$mu, lwd = 3, type= "l", col = "blue")
  # lines(t_grid, fit.lspensp$mu, lwd = 3, type = "l", col = "red")
  
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


################################################################################################

library(microbenchmark)

## --- Simulation with heavy-tailed KL process --- ##
nsim <- 200
n    <- 100
p    <- 50

mse.smsp <- mse.pensp <- mse.lspensp <- rep(0, nsim)
time.smsp <- time.pensp <- time.lspensp <- rep(0, nsim)

shapes.smsp <- shapes.pensp <- shapes.lspensp <- matrix(0, ncol = nsim, nrow = p)

t_grid <- seq(0, 1, length.out = p)
Y <- matrix(NA, nrow = n, ncol = p)

## Population mean
# mu_true <- function(t) exp(-(t-0.25)^2/0.01) +
#   exp(-(t-0.50)^2/0.01) +
#   exp(-(t-0.75)^2/0.01) 
mu_true <- function(t) sin(2*pi*t)
mu_grid <- mu_true(t_grid)

# set.seed(123)

for(k in 1:nsim){
  print(k)
  
  ## --- Generate data ---
  X <- matrix(NA, nrow = n, ncol = p)
  for(i in 1:n){
    X[i,] <- mu_grid 
    for(j in 1:50){ 
      X[i,] <- X[i, ] + sqrt(2)*rt(1, df = 5) *
        sapply(t_grid, function(x) sin((j-1/2)*pi*x)/((j-1/2)*pi) )
    }
    m_i <- sample(floor(0.5 * p):floor(0.8 * p), 1)
    idx <- sort(sample(seq_len(p), m_i))
    zeta <- 0.5*rt(m_i, df = 2)
    Y[i, idx] <- X[i, idx] + zeta
  }
  
  ## --- Fit estimators and record times ---
  # t1 <- system.time(fit.smsp <- quan_smsp(Y, alpha = 0.5))
  # t2 <- system.time(fit.pensp <- quan_pensp(Y, alpha = 0.5, K = 30))
  # t3 <- system.time(fit.lspensp <- ls_pensp(Y, K = 30))
  # 
  # time.smsp[k]   <- t1["elapsed"]
  # time.pensp[k]  <- t2["elapsed"]
  # time.lspensp[k] <- t3["elapsed"]
  t1 <- microbenchmark(fit.smsp <- quan_smsp(Y, alpha = 0.5), times = 3)
  t2 <- microbenchmark(fit.pensp <- quan_pensp(Y, alpha = 0.5), times = 3)
  # t2 <- microbenchmark(fit.lspensp <- quan_smsp(Y, alpha = 0.5), times = 3)
  
  time.smsp[k]   <- mean(t1$time/1e+09)
  time.pensp[k]  <- mean(t2$time/1e+09)
  
  ## --- Compute MSEs ---
  mse.smsp[k]   <- mean((fit.smsp$mu - mu_grid)^2)
  mse.pensp[k]  <- mean((fit.pensp$mu - mu_grid)^2)
  # mse.lspensp[k] <- mean((fit.lspensp$mu - mu_grid)^2)
  
  ## --- Save fitted shapes ---
  shapes.smsp[, k]   <- fit.smsp$mu
  shapes.pensp[, k]  <- fit.pensp$mu
  # shapes.lspensp[, k] <- fit.lspensp$mu
}

## --- Summarise results ---
results <- data.frame(
  Estimator = c("SmSP", "PenSP", "LS-PenSP"),
  MSE_mean  = 1000*c(mean(mse.smsp), mean(mse.pensp), mean(mse.lspensp)),
  MSE_se    = 1000*c(sd(mse.smsp)/sqrt(nsim),
                sd(mse.pensp)/sqrt(nsim),
                sd(mse.lspensp)/sqrt(nsim)),
  Time_mean = c(mean(time.smsp), mean(time.pensp), mean(time.lspensp)),
  Time_se   = c(sd(time.smsp)/sqrt(nsim),
                sd(time.pensp)/sqrt(nsim),
                sd(time.lspensp)/sqrt(nsim))
)

print(results)
## --- Optional: boxplots of times and MSEs ---
par(mfrow=c(1,2))
boxplot(time.smsp, time.pensp, time.lspensp,
        names=c("SmSP", "PenSP", "LS-PenSP"),
        main="Runtime comparison", ylab="Seconds")
boxplot(mse.smsp, mse.pensp, mse.lspensp,
        names=c("SmSP", "PenSP", "LS-PenSP"),
        main="MSE comparison", ylab="MSE")

## --- Optional: boxplots of times and MSEs ---
par(mfrow=c(1,1))
# save(time.smsp, time.pensp, file = "df2.RData")
par(mar = c(4,5.5,2,2), mgp = c(3.6, 1.5, 0))
boxplot(log(time.smsp), log(time.pensp), pch = 20, cex = 2, lwd = 1.5,
        names=c("SmSP", "PenSP"), ylab="log(Seconds)", cex.axis = 2.5, cex.lab = 2.5);grid()

###############################################################

files <- c("df1.RData", "df2.RData", "df5.RData", "df10.RData")
df_labels <- c("df1", "df2", "df5", "df10")

Estimator <- character()
Time <- numeric()
DF <- character()

for (i in seq_along(files)) {
  vars <- load(files[i])  # loads time.smsp, time.pensp, time.lspensp
  
  # Append results
  Estimator <- c(Estimator,
                 rep("SmSP", length(time.smsp)),
                 rep("PenSP", length(time.pensp)))
  
  Time <- c(Time, time.smsp, time.pensp)
  
  DF <- c(DF,
          rep(df_labels[i], length(time.smsp) + length(time.pensp)))
  
  rm(list = vars)
}

all_times <- data.frame(Estimator, Time, DF)

# Side-by-side grouped boxplot
par(mfrow=c(1,1))
par(mar = c(5,5.5,2,2), mgp = c(3.6, 1.5, 0))
boxplot(log(Time) ~ Estimator + DF, data = all_times, lwd = 1.5, cex.axis = 2, cex.lab = 2.5, pch = 20,
        las = 2, xlab = "", xaxt = "n", # rotate labels
        col = rep(c("skyblue", "salmon"), times = length(df_labels)),
        ylab = "log(seconds)")
grid()

# Small estimator labels under each box
est_labels <- rep(c("PenLAD", "SmLAD"), times = length(df_labels))
axis(1, at = 1:length(est_labels), labels = est_labels, tick = FALSE, line = 0.1, cex.axis = 1.8)

# Big DF labels centered under each pair of boxes
df_positions <- sapply(1:length(df_labels), function(i) mean((2*i-1):(2*i)))  # centers
axis(1, at = df_positions, labels = c("df = 1", "df = 2", "df = 5", "df = 10"), tick = FALSE, line = 2, 
     cex.axis = 2.2)
