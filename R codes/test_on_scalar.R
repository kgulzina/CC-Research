####### This code runs GP model with scalar input and output #######
## log is given "README.md" file in github
## Author: stat_cat
## Date: 05/06/19 - present




# Libraries ---------------------------------------------------------------
library(MASS)
library(mvtnorm)
library(Matrix)




# Conventions -------------------------------------------------------------
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





# Model -------------------------------------------------------------------
### write here the loglkl function 



### write here the omega function



### simulate w






# Estimation --------------------------------------------------------------































