####### This code runs GP model with scalar input and output #######
## log is given "README.md" file in github
## Author: stat_cat
## Date: 05/06/19 - present




# Libraries ---------------------------------------------------------------
library(MASS)
library(mvtnorm)
library(Matrix)




# Conventions -------------------------------------------------------------
n <- 10 # for more sample size I am getting more stable estimates!
rows <- sample(1:432, n)
q <- 1
rho <- 0.99
sigmasq <- 1

## to use supplementary functions:
source('R codes/supplementary_functions.R')
#source('R codes/loglkl_with_penalty.R')
#source("R codes/wepp_data.R")




# Data --------------------------------------------------------------------
d <- read.csv("data/final_csv/annual_soil_loss_not_stdrzd.csv")
d_sf <- data.frame(total_prcp = apply(d[,5:369], 1, sum),
                       d[,370:384],
                       soil_loss = d[,4])

# standardize precipitation amount
d_sf[,1] <- log(d_sf[,1] + min(d_sf[d_sf[,1] != 0,1]))



# Model -------------------------------------------------------------------
### write here the omega function
gcalc_corr_scalar <- function(d,w) {
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
    y <- ncol(d)
    omega <- matrix(NA, nrow = n, ncol = n)
    
    for(i in 1:n){
        for(j in 1:n){
            omega[i,j] <- exp(sum(-w*(d[i,-y] - d[j,-y])^2))
        }
    }
    return(omega)
}



### write here the loglkl function 
lkl_sf_gp <- function(theta, d) {
    #
    # Args:
    #   theta:  parameters = (w1, w2, lambda) in a vector form
    #   d:      data frame, last column is response Y, others are input X's
    #
    # Output:
    #   result: likelihood value
    
    # constants
    n <- nrow(d)
    thetal <- length(theta)
    q <- (thetal-2) / 2
    
    # penalty on w
    basis <- generate_trig_basis(14, q)
    logit_w <- theta[2:(thetal-1)]%*%basis
    
    # get w
    w <- exp(logit_w)/(1+exp(logit_w))
    
    # calculate the covariance matrix for yw
    omega <- gcalc_corr_scalar(d, c(theta[1], w)) + diag(x = theta[thetal],
                                                     nrow = n,
                                                     ncol = n)
    
    # data model: log_likelihood 
    p_y <- dmvnorm(d[,ncol(d)], mean = rep(0, n), sigma = omega)  
    
    # the posterior which will be maximized
    result <- p_y
    
    return(result)
}



# Estimation --------------------------------------------------------------
pars_theta <- c(0.5, 0.5, 0.2, 0.1, 0.1)
sf_gp_est <- optim(par = pars_theta, 
                       lkl_sf_gp, 
                       d = d_sf[rows,], 
                       control = list(fnscale = -1,
                                      maxit = 1000))
sf_gp_est



# results
# n = 20, pars_theta <- c(0.5, 0.5, 0.2, 0.1, 0.1)
# >>> par = 2.174854e-02 -5.440885e+02  1.069658e+03 -2.776159e+02  
# 4.637701e-03) CVG = 0

# n = 50, pars_theta <- c(0.5, 0.5, 0.2, 0.1, 0.1)
# >>> par = (0.593970780 -0.623168211  0.003892142) , CVG = 10


# n = 100, pars_theta <- c(0.5, 0.5, 0.2, 0.1, 0.1)
# >>> par = (0.447994602 -0.095504035  0.003866797), CVG = 10

# n = 150, pars_theta <- c(0.5, 0.5, 0.2, 0.1, 0.1)
# >>> par = (0.545756609 -0.018987352  0.003842057)

# n = 200, pars_theta <- c(0.5, 0.5, 0.2, 0.1, 0.1)
# >>> par = (0.555937239 -0.007558199  0.003858029), CVG = 0

# n = 300, pars_theta <- c(0.5, 0.5, 0.2, 0.1, 0.1)
# >>> par = ()

# n = 432, pars_theta <- c(0.5, 0.5, 0.2, 0.1, 0.1)
# >>> par = ()