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
n <- 120 # for more sample size I am getting more stable estimates!
Time <- 364
q <- 1 #number of harmonics
rho <- 0.99
sigmasq <- 1

## to use supplementary functions:
source('supplementary_functions.R')
source('loglkl_with_penalty.R')
#source("wepp_data.R")




# Data --------------------------------------------------------------------
d <- read.csv("data/annual_soil_loss.csv") #home/kgulzina/cc/
d <- d[1:n, c(2:366, 369)]

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