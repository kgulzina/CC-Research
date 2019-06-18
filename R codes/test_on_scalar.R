####### This code runs GP model with scalar input and output #######
## log is given "README.md" file in github
## Author: stat_cat
## Date: 05/06/19 - present




# Libraries ---------------------------------------------------------------
library(MASS)
library(mvtnorm)
library(Matrix)




# Conventions -------------------------------------------------------------
n <- 150 # for more sample size I am getting more stable estimates!
rows <- sample(1:432, n)
rho <- 0.99
sigmasq <- 1

## to use supplementary functions:
source('R codes/supplementary_functions.R')
source('R codes/loglkl_with_penalty.R')
source("R codes/wepp_data.R")




# Data --------------------------------------------------------------------
d <- read.csv("data/final_csv/annual_soil_loss_not_stdrzd.csv")
d_scalar <- data.frame(total_prcp = apply(d[,5:369], 1, sum),
                       avg_slope = apply(d[370:384], 1, mean),
                       soil_loss = d[,4])

# standardize precipitation amount
d_scalar[,1] <- log(d_scalar[,1] + min(d_scalar[d_scalar[,1] != 0,1]))
d_scalar[,2] <- d_scalar[,2] / 100



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
lkl_scalar_gp <- function(theta, d) {
    #
    # Args:
    #   theta:  parameters = (w1, w2, lambda) in a vector form
    #   d:      data frame, last column is response Y, others are input X's
    #
    # Output:
    #   result: likelihood value
    
    # constants
    n <- nrow(d)
    
    # calculate the covariance matrix for yw
    omega <- gcalc_corr_scalar(d, theta[1:2]) + diag(x = theta[3],
                                           nrow = n,
                                           ncol = n)
    
    # data model: log_likelihood 
    p_y <- dmvnorm(d[,ncol(d)], mean = rep(0, n), sigma = omega)  
    
    # the posterior which will be maximized
    result <- p_y
    
    return(result)
}



# Estimation --------------------------------------------------------------
pars_theta <- c(0.5 , 0.5, 0.1)
scalar_gp_est <- optim(par = pars_theta, 
                       lkl_scalar_gp, 
                       d = d_scalar[1:300,], 
                       control = list(fnscale = -1,
                                      maxit = 1000))
scalar_gp_est



# results
# n = 20, pars_theta <- c(0.5 , 0.5, 0.1) 
# >>> par = (0.080920239 -1.142949670  0.004392549), CVG = 0

# n = 50, pars_theta <- c(0.080920239 -1.142949670  0.004392549)
# >>> par = (0.593970780 -0.623168211  0.003892142) , CVG = 10


# n = 100, pars_theta <- c(0.5 , 0.5, 0.1)
# >>> par = (0.447994602 -0.095504035  0.003866797), CVG = 10

# n = 150, pars_theta <- c(0.5 , 0.5, 0.1)
# >>> par = (0.545756609 -0.018987352  0.003842057)

# n = 200, pars_theta <- c(0.5 , 0.5, 0.1)
# >>> par = (0.555937239 -0.007558199  0.003858029), CVG = 0

# n = 300, pars_theta <- c(0.5 , 0.5, 0.1)
# >>> par = ()

# n = 432, pars_theta <- c(0.5 , 0.5, 0.1)
# >>> par = ()































