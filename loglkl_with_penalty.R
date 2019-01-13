####### Log-likelihood function for GP with penalty on weight: w #######
library(MASS)
library(mvtnorm)
library(Matrix)
library(tmvtnorm)






# cov. matrix for MVN simulation: n samples
#sigma <- calc_Sigma(Time+1, rho = 0.95)

# simulate X
#x_values <- mvrnorm(n, rep(0, Time+1), sigma) # vector of length T+1

# a correlation matrix: deterministic function for w(t)
#calc_corr <- function(val, t, n) { 
    # Define w(t)=(T-t)/T :straight line with the negative slope: -1/T
    # Calculate covariances between simulated X's using correlation form 
    # on page3, (1) form.
    #omega <- matrix(NA, nrow = n, ncol = n)
    #time <- 0:t
    #for(i in 1:n) {
     #   for(j in 1:n) {
      #      omega[i,j] <- exp(-sum(((t-time+1)/(t+1))*(val[i,]-val[j,])^2))
       # }
    #}
    #return(omega)
#}

#omega <- calc_corr(x_values, Time, n)


# simulate W's
#simulate_w <- function(t){
 #   sigma <- 1/(0.36)*calc_Sigma(t+1, 0.99)
    # simulate from truncated normal
  #  w <- rtmvnorm(1, rep(0, t+1), sigma, lower = rep(0, t+1), 
   #               upper = rep(Inf, t+1), algorithm = "gibbs")
    #w <- mvrnorm(1, rep(0, t+1), sigma)
    #return(w)
#}

#w <- simulate_w(Time)

# a correalation matrix: w(t) sampled from AR(1)
#calc_corr2 <- function(val, t, n, w){
    
    # Weights are sampled
 #   omega <- matrix(NA, nrow = n, ncol = n)
  #  time <- 0:t
   # for(i in 1:n) {
    #    for(j in 1:n) {
     #       omega[i,j] <- exp(-sum(w*(val[i,]-val[j,])^2))
      #  }
    #}
    #return(omega)
    
    
    
#}

#omega2 <- calc_corr2(x_values, Time, n, w)


# simulate y
#y <- rmvnorm(1, mean = rep(0, n), sigma = omega2) 



# d: dataset combined: functional input and scale output
#d <- cbind(x_values,t(y)) #colnames(d) <- c("t0", "t1", "y")



# need omega -- general w/ w(t) as a vector. 
#gcalc_corr <- function(d,w) { #retriev x's from d
 #   n <- nrow(d)
  #  time <- ncol(d)-1
   # omega <- matrix(NA, nrow = n, ncol = n)
    
    # remove last column: y
    #for(i in 1:n){
     #   for(j in 1:n){
      #      omega[i,j] <- exp(-sum(w*(d[i,-ncol(d)]-d[j,-ncol(d)])^2))
       # }
    #}
    #return(omega)
#}

#pars <- seq(1, 1/(Time+1), len = Time+1)
#trial <- gcalc_corr(d, pars) 
#rankMatrix(trial)



# new log_lkl w/ penalty term to optimize:
loglkl_mvn_penalty <- function(w,d) { # 
    
    time <- ncol(d)-2
    # calculate the covariance matrix for yw
    omega <- gcalc_corr(d,w)
    
    # calculate the convariance matrix for w
    sigma <- 1/(0.36)*calc_Sigma(time+1, 0.99)
    
    # data model: log_likelihood 
    p_yw <- dmvnorm(d[,ncol(d)], mean = rep(0, nrow(d)), sigma = omega, log = TRUE)  
    
    # prior on w: log_likelihood
    p_w <- dmvnorm(w, mean = rep(0, time+1), sigma = sigma , log = TRUE)
    
    # the posterior which will be maximized
    result <- p_yw + p_w 
    
    return(result)
}





gradient_loglkl_penalty <- function(w, d){
# Calculates the gradient of f: log-likelihood with penalty w.r.t w
#
# Args: 
#   w: weights, in a vector form
#   d: data frame, last column is response Y, others are input X's
#
# Output:
#   gr: gradients, in a vector form
    
    n <- length(w)
    gr <- numeric(n)
    
    # omega
    omega <- calc_corr_sampled(d[,-(n+1)], n-1, nrow(d), w)
    
    # inverse of omega
    inv_omega <- solve(omega)
    
    # inverse of pi
    inv_pi <- solve(1/(0.36)*calc_Sigma(n, 0.99))
    
    for(i in 1:n) {
        dm <- calc_mtrx_deriv(d[,-ncol(d)], n-1, nrow(d), w, i)
        # for penalty term
        dpenalty <- 0
        for (j in 1:n){
            temp = w[j]*(inv_pi[i,j] + inv_pi[j,i])
            dpenalty <- dpenalty + temp
        }
        
        gr[i] <- -1/2*sum(diag(inv_omega%*%dm)) + 
                  1/2*t(d[,n+1])%*%inv_omega%*%dm%*%inv_omega%*%d[,n+1] +
                  -1/2*dpenalty
    }
    
    return(gr)
                                                                                                                                                                                               
}



calc_mtrx_deriv <- function(X, t, n, w, k){
# Calculates derivatives of matrix w.r.t one specific parameter
#
# Args:
#   X: simulated, in a vector form
#   t: time(T) or length(X) - 1
#   n: sample size
#   w: simulated weights, in a vector form 
#   k: index of parameter vector, by which derivative is taken
#
# Output:
#   dm: matrix of elementwise derivatives
    
    dm <- matrix(NA, nrow = n, ncol = n)
    time <- 0:t
    for(i in 1:n) {
        for(j in 1:n) {
            dm[i,j] <- -exp(-sum(w*(X[i,]-X[j,])^2))*((X[i,k]-X[j,k])^2)
        }
    }
    return(dm)
}



# use log_likelihood w/ penalty:
#opt <- optim(par = pars, loglkl_mvn_penalty, d = d,
 #             control = list(fnscale = -1, maxit=4000))
#opt$convergence 
#opt$par
#pars


# compare w/ pars2:
#opt2 <- optim(par = pars2, loglkl_mvn_penalty, d = d,
 #            control = list(fnscale = -1, maxit=4000))
#opt2$convergence 
#opt2$par
#pars2




######## Automation ##########

# constant terms: simuated data
#simulate_d <- function(Time, n, w){
    
 #   sigma <- Sigma(Time+1, rho = 0.95)
    
    # simulate X's
  #  x_values <- mvrnorm(n, rep(0, Time+1), sigma)
    
    #calculate omega for Y
   # omega <- calc_corr2(x_values, Time, n, w)
    
    # simulate
    #y <- rmvnorm(1, mean = rep(0, n), sigma = omega)
    
    # combine datasets: functional input and scale output
    #d <- cbind(x_values,t(y))
    
    #return(d)
    
#}

#w <- simulate_w(10) # time 
#d <- simulate_d(10, 100, w) # time, n


# constant term: pars
#pars <- seq(1, 1/(10+1), len = 10+1) # for input


#estimate_w <- function(opt_f, grad, pars, d) { 
    
    # returns optim() output with tru parameter values
 #   opt <- optim(par = pars, opt_f, d = d, control = list(fnscale = -1,
  #                                                   maxit=10000),
   #              gr = grad) # added gradient function
    
#    print(pars)
 #   print(w)
  #  return(opt)
#}



# comments:
#### new concerns.. should get rid of negative weights?
#### truncate or not truncated? >>> have to find new density functions
#### trying 365 on server:  Error in checkSymmetricPositiveDefinite(sigma) : 
#sigma must be positive definite ??? 


































