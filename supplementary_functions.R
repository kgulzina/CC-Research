####### Complementary/supplementary fuctions #######
## Author: stat_cat
## Date: 01/03/19

calc_Sigma <- function(n, rho){ 
# Calculates the covariance matrix of AR(1) process
#  
# Args:
#   n:   vector(sample) size (integer)
#   rho: fixed correlation value (|real|<=1)
#
# Output:
#   Sigma: symmetric covariance matrix
    
    Sigma <- matrix(NA, nrow = n, ncol = n)
    for(i in 1:n){
        for(j in 1:n){
            Sigma[i,j] <- rho^abs(i-j)
        }
    }
    return(Sigma)
}




calc_corr_deterministic <- function(X, t, n) { 
# Calculate covariances between simulated X's using correlation form used in 
# the GP model assumptions: deterministic w(t)
#
# Args: 
#   X: simulated, in a vector form
#   t: time(T) or length(X) - 1
#   n: sample size
#
# Output:
#   omega: correlation matrix of size n 
    
    omega <- matrix(NA, nrow = n, ncol = n)
    time <- 0:t
    for(i in 1:n) {
        for(j in 1:n) {
            omega[i,j] <- exp(-sum(((t-time+1)/(t+1))*(X[i,]-X[j,])^2))
        }
    }
    return(omega)
}




calc_corr_sampled <- function(X, t, n, w){
# Calculates covariances between simulated X's using correlation form used in 
# the GP model assumptions: sampled w(t) from specific density
#
# Args: 
#   X: simulated, in a vector form
#   t: time(T) or length(X) - 1
#   n: sample size
#   w: simulated weights, in a vector form
#
# Output:
#   omega: correlation matrix of size n    
    
    omega <- matrix(NA, nrow = n, ncol = n)
    time <- 0:t
    for(i in 1:n) {
        for(j in 1:n) {
            omega[i,j] <- exp(-sum(w*(val[i,]-val[j,])^2))
        }
    }
    return(omega)
}




gcalc_corr <- function(d,w) { #retriev x's from d
# Calculates covariances between X's using correlation form used in the GP
# model assumptions: to find MLE of w(t)
#
# Args: 
#   d: data frame, last column is response Y, others are input X's
#   w: weights, in a vector form
#
# Output:
#   omega: correlation matrix of size n
    
    n <- nrow(d)
    time <- ncol(d)-1
    omega <- matrix(NA, nrow = n, ncol = n)
    
    for(i in 1:n){
        for(j in 1:n){
            omega[i,j] <- exp(-sum(w*(d[i,-ncol(d)]-d[j,-ncol(d)])^2))
        }
    }
    return(omega)
}




sim_truncated_w <- function(t){
# Simulated weights in covariance function of GP (not deterministic) from 
# truncated MVN (see assumptions in "CC-Process" log-file)
#
# Args: 
#   t: time(T) or length(X)-1
#
# Output:
#   w: weights, in a vector form
    
    sigma <- 1/(0.36)*Sigma(t+1, 0.99)
    # simulate from specific density
    w <- rtmvnorm(1, rep(0, t+1), sigma, lower = rep(0, t+1), 
                  upper = rep(Inf, t+1), algorithm = "gibbs")
    # w <- mvrnorm(1, rep(0, t+1), sigma)
    return(w)
}




calc_gradient_num <- function(f,w,d,epsilon=10^-8){
# Calculates the gradient of a function numerically
# 
# Args:
#   f: function (log-likelihood function)
#   d: data frame, last column is response Y, others are input X's
#   w: weights, in a vector form 
#   
# Output:
#   gr: gradients, in a vector form
    
    n <- length(w)
    gr <- numeric(n) 
    for(i in 1:n) {
        h <- rep(0,n); h[i] <- epsilon
        gr[i] <- (f(w+h,d)-
                      f(w,d))/epsilon
    }
    return(gr)
}




simulate_d <- function(t, n, w){
# Simulates data: both X and Y according to GP model described in 
# "CC-Process" log-file
# 
# Args:
#   t: Time(T) or length(X)-1
#   n: sample size
#   w: weights, in a vector form
#
# Outputs:
# d: data frame, last column is response Y, others are input X's 
    
    sigma <- Sigma(t+1, rho = 0.95)
    
    # simulate X's
    x_values <- mvrnorm(n, rep(0, t+1), sigma)
    
    #calculate omega for Y
    omega <- calc_corr2(x_values, t, n, w)
    
    # simulate Y's
    y <- rmvnorm(1, mean = rep(0, n), sigma = omega)
    
    # combine datasets: functional input and scale output
    d <- cbind(x_values,t(y))
    
    return(d)
}




simulate_w <- function(t){
# 
#
#
#
#
    sigma <- 1/(0.36)*Sigma(t+1, 0.99)
    # simulate from truncated normal
    w <- rtmvnorm(1, rep(0, t+1), sigma, lower = rep(0, t+1), 
                  upper = rep(Inf, t+1), algorithm = "gibbs")
    #w <- mvrnorm(1, rep(0, t+1), sigma)
    return(w)
}







































