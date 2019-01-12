######### Complementary/supplementary fuctions #########
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
            omega[i,j] <- exp(-sum(w*(X[i,]-X[j,])^2))
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
    
    sigma <- calc_Sigma(t+1, rho = 0.95)
    
    # simulate X's
    x_values <- mvrnorm(n, rep(0, t+1), sigma)
    
    #calculate omega for Y
    omega <- calc_corr_sampled(x_values, t, n, w)
    
    # simulate Y's
    y <- rmvnorm(1, mean = rep(0, n), sigma = omega)
    
    # combine datasets: functional input and scale output
    d <- cbind(x_values,t(y))
    
    return(d)
}




simulate_w_mvn <- function(t){
# Simulates weights from MVN ~ (GP), according to assumption in 
# CC-Process - 1c
#
# Args:
#   t: Time(T) or length(X)-1
#
# Output:
#   w: weights, in a vector form
    
    sigma <- 1/(0.36)*calc_Sigma(t+1, 0.99)
    # sigma depends on t, so not constant
    w <- mvrnorm(1, rep(0, t+1), sigma)
    return(w)
}




simulate_w_trnctd <- function(t){
# Simulates weights from truncated (positive) MVN ~ (GP), according to 
# the assumptions in CC-Process - 1c
#    
# Args:
#   t: Time(T) or length(X)-1
#
# Output:
#   w: weights, in a vector form
    
    sigma <- 1/(0.36)*calc_Sigma(t+1, 0.99)
    # sigma depends on t, so not constant
    w <- rtmvnorm(1, rep(0, t+1), sigma, lower = rep(0, t+1), 
                  upper = rep(Inf, t+1), algorithm = "gibbs")
    return(w)
}




simulate_log_w_stan <- function(t){ ## this function does not work right know
# Simulates log(weights) from transformed MVN ~ (GP) given in CC-Process 1c  
# and returns transformed result
#    
# Args: 
#   t: Time(T) or length(X)-1
#
# Output:
#   w: weights, in a vector form
    
    # define the density
    sigma <- 1/(0.36)*calc_Sigma(t+1, 0.99)
    dens <- function(theta){
        # theta is a vector of length T
        result <- (2*pi)^(-(t+1)/2)*(det(sigma))^(-1/2)*exp(sum(theta))*
            exp(-1/2*theta%*%solve(sigma)%*%t(theta))
        return(result) ### need a scalar, not vector
    }
    
    require(mcmc)
    out <- metrop(dens, sigma, 1e+3)
    return(out$batch)
}




simulate_log_w <- function(t){
# Simulates log(weights) from  MVN(-2, Pi) given in CC-Process 1c and  
# returns transformed result
#    
# Args: 
#   t: Time(T) or length(X)-1
#
# Output:
#   w: weights, in a vector form
    
    sigma <- 1/(0.36)*calc_Sigma(t+1, 0.99)
    # sigma depends on t, so not constant
    theta <- mvrnorm(1, rep(-3, t+1), sigma)
    w <- exp(theta)
    return(w)
}





estimate_w <- function(opt_f, grad, pars, d, maxit){ 
# Estimates w using eBayes approach, i.e finds MLE estimates of w.
#    
# Args:   
#   opt_f: likelihood function to be optimized
#   grad:  gradient of opt_f
#   pars:  initial values for w
#   d:     observed data
#   maxit: maximum number of iterations
#
# Output: 
#   opt:   results of optim() 
    
    opt <- optim(par = pars, opt_f, d = d, control = list(fnscale = -1,
                                                          maxit=maxit),
                 gr = grad) 
    
    # print true values of w
    print("True values of w")
    print(w)
    
    # estimated w
    return(opt)
}







































