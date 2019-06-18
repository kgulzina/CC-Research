####### This code tests "loglkl_with_penalty" with the actual data #######
## log is given "README.md" file in github
## Author: stat_cat
## Date: 09/03/18 - present



# Libraries ---------------------------------------------------------------
library(MASS)
library(mvtnorm)
library(Matrix)
#library(dlm)




# Conventions -------------------------------------------------------------
## sample size and Time: T << n
#### warning!!! Time = actual time - 1 = 365 - 1 = 364
#n <- 120 # for more sample size I am getting more stable estimates!
#Time <- 364
q <- 1 #number of harmonics
rho <- 0.99
sigmasq <- 1

## to use supplementary functions:
source('R codes/supplementary_functions.R')
source('R codes/loglkl_with_penalty.R')
#source("wepp_data.R")




# Data --------------------------------------------------------------------
d <- read.csv("data/final_csv/annual_soil_loss.csv") #home/kgulzina/cc/
d <- cbind(d[,-c(1:4)], d[, 4])

### commments: here first 365 columns are for precipitation in (mm),
# the las columnt is output = annual soil loss in (kg/m^2).




# Initial values ----------------------------------------------------------
### should work on these initial values!!!
pars_coef <- simulate_coeff_dlm(q)
pars_coef <- rep(0, 2*q+1)
pars_coef <- rep(1, 2*q+1)
pars_coef
# try another initial values from truncated n(0,1)
# try another initial values from gamma()
# try another initial values from beta()

# to know what is going on:
dynamic_loglkl_mvn_penalty(pars_coef, d)
### Comments: since the climate file is the same for all 12*3 hills, we have 
### all etnries of omega == 1, when calculating correlation matrix!!!


# Dynamic Linear Models ---------------------------------------------------
est3_gr1_pars <- optim(par = pars_coef, 
                    dynamic_loglkl_mvn_penalty, 
                    d = d, 
                    control = list(fnscale = -1,
                    maxit = 1000000))
est3_gr1_pars
theta




# Two numerical inputs ----------------------------------------------------
## we need new dynamic_loglkl_mvn_penalty function
dynamic_lkl_mvn_penalty2 <- function(theta, d) {
    #
    # Args:
    #   theta:
    #   d:
    #
    # Output:
    #   result:
    
    # some constants
    time1 <- ncol(d)-17
    time2 <- 14 #fixed for now
    
    # length of theta parameter
    thetalength <- length(theta) / 2
    
    # divide coefficients into two
    theta1 <- theta[1:thetalength]
    theta2 <- theta[-(1:thetalength)]
     
    # q is fixed for both   
    q <- (length(theta1) - 1) / 2
    
    # penalty on w
    basis1 <- generate_trig_basis(time1+1, q)
    logit_w <- theta1%*%basis1
    
    # get w
    w <- exp(logit_w)/(1+exp(logit_w))
    
    # penalty on v
    basis2 <- generate_trig_basis(time2+1, q)
    logit_v <- theta2%*%basis2
    
    # get v
    v <- exp(logit_v)/(1+exp(logit_v))
    
    # combine two weights into one
    wv <- c(w,v)
    
    # calculate the covariance matrix for ywv
    omega <- gcalc_corr(d,wv)
    
    # data model: log_likelihood 
    p_ywv <- dmvnorm(d[,ncol(d)], mean = rep(0, nrow(d)), sigma = omega)  
    
    # the posterior which will be maximized
    result <- p_ywv
    
    return(result)
}



# Initial values for double GP input model ---------------------------------
theta <- c(0.05, 0.002, 0.002)
pars_coef <- c(theta, theta)
pars_coef



# Run double input GP model -----------------------------------------------
n <- 20
rows <- sample(1:432, n)
fun2_gp_est <- optim(par = pars_coef, 
                       dynamic_lkl_mvn_penalty2, 
                       d = d[rows,], 
                       control = list(fnscale = -1,
                                      maxit = 100))
fun2_gp_est
theta







# One scalar, one functional input ----------------------------------------
## we need new dynamic_loglkl_mvn_penalty function
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

# results
# n = 20, pars_theta <- c(0.5 , 0.5, 0.1) 
n <- 20
rows <- sample(1:432, n)
sf_gp_est <- optim(par = pars_theta, 
                   lkl_sf_gp, 
                   d = d_sf[rows,], 
                   control = list(fnscale = -1,
                                  maxit = 1000))
sf_gp_est








# Test on wepp data with penalty ------------------------------------------
dynamic_lkl_mvn_penalty_ridge <- function(theta, d) {
    #
    # Args:
    #   theta:
    #   d:
    #
    # Output:
    #   result:
    
    # some constants
    n <- nrow(d)
    time1 <- ncol(d)-17
    time2 <- 14 #fixed for now
    
    # length of theta parameter
    lambdal <- length(theta)
    thetalength <- (lambdal - 1) / 2
    
    
    # divide coefficients into two
    theta1 <- theta[1:thetalength]
    theta2 <- theta[-c(1:thetalength,lambdal)]
    
    # q is fixed for both   
    q <- (length(theta1) - 1) / 2
    
    # penalty on w
    basis1 <- generate_trig_basis(time1+1, q)
    logit_w <- theta1%*%basis1
    
    # get w
    w <- exp(logit_w)/(1+exp(logit_w))
    
    # penalty on v
    basis2 <- generate_trig_basis(time2+1, q)
    logit_v <- theta2%*%basis2
    
    # get v
    v <- exp(logit_v)/(1+exp(logit_v))
    
    # combine two weights into one
    wv <- c(w,v)
    
    # calculate the covariance matrix for ywv + penalty
    omega <- gcalc_corr(d,wv) + diag(theta[lambdal],
                                     ncol = n,
                                     nrow = n)
    
    # data model: log_likelihood 
    p_ywv <- dmvnorm(d[,ncol(d)], mean = rep(0, nrow(d)), sigma = omega)  
    
    # the posterior which will be maximized
    result <- p_ywv
    
    return(result)
}



# Initial values for double GP input model ---------------------------------
pars_coef <- c(0.5, 0.1, 0.1, 0.5, 0.1, 0.1, 5)
#pars_coef <- c(simulate_coeff_dlm(q), simulate_coeff_dlm(q))
pars_coef



# Run double input GP model -----------------------------------------------
n <- 20
rows <- sample(1:432, n)
fun2_gp_est <- optim(par = pars_coef, 
                     dynamic_lkl_mvn_penalty_ridge, 
                     d = d[rows,], 
                     control = list(fnscale = -1,
                                    maxit = 1000))
fun2_gp_est














