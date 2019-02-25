####### This code tests "loglkl_with_penalty" with the simulated data #######
## log is given "README.md" file in github
## Author: stat_cat
## Date: 09/03/18 - present




# Libraries ---------------------------------------------------------------
library(MASS)
library(mvtnorm)
library(Matrix)
library(tmvtnorm)
library(plotly)
library(lattice)
library(dlm)




# Conventions -------------------------------------------------------------
## sample size and Time: T << n
n <- 12 # for more sample size I am getting more stable estimates!
Time <- 10
rho <- 0.99
sigmasq <- 1

## to use supplementary functions:
source('supplementary_functions.R')
source('loglkl_with_penalty.R')
source("assess_accuracy.R")
source("3D_plot_loglkl.R")




# Simulation --------------------------------------------------------------
## simulate weights using truncated mvn:
#w <- simulate_w_trnctd(Time)

## Approach 1 to simulate w:
#w1 <- simulate_log_w_stan(Time)

## Approach 2 to simulate w:
w <- simulate_log_w(Time) #current

## simulate data:
d <- simulate_d(Time, n, w)

## simulate truncated data:
d <- simulate_trunc_d(Time, n, w)

## real wepp data
d <- wepp_data[,c(1:11, ncol(wepp_data))]



# If Y’s are sampled correctly --------------------------------------------
## summary statistics and simple plot
y <- d[,ncol(d)]
summary(y)
plot(y)

## color map:
# we need positive y_values:
positive_y <- y + abs(min(y)) #shift by the minimum value
plot(positive_y)
min1_y <- positive_y + (1-min(positive_y)) #set minimum value of y = 1
plot(min1_y)

# convert to colors:
colors <- colorRampPalette(c("red", "orange", "blue"))(3)[min1_y]
 

# check if all Y's have been assigned colors: n is good
length(colors)

# add to dataframe:
to_check <- as.data.frame(d[,c(1,2)])
colnames(to_check) <- c("x0","x1")
to_check$colors <- colors

# plotted are X's (default: n=100) and points are colored by Y:
plot(to_check$x0, to_check$x1, col = to_check$colors)

##### close dots (dy distance) should be of the same color #####




# Initial values ----------------------------------------------------------
## set initial values as strict line:
pars1 <- seq(1, 1/(Time+1), len = Time+1)
pars1

## set initial values sampled from AR(1):
pars2 <- simulate_log_w(Time)
pars2




# Estimation --------------------------------------------------------------
## compare two estimates with different initial values using num gradient:
est_gr1_pars1 <- estimate_w(loglkl_mvn_penalty, calc_gradient_num, pars1, d, 100000)
est_gr1_pars1

est_gr1_pars2 <- estimate_w(loglkl_mvn_penalty, calc_gradient_num, pars2, d, 100000)
est_gr1_pars2

est_gr1_w <- estimate_w(loglkl_mvn_penalty, calc_gradient_num, w, d, 100000)
est_gr1_w #initial values as true parameters


## compare two estimates with different initial values and true gradient:
est_gr2_pars1 <- estimate_w(loglkl_mvn_penalty, gradient_loglkl_penalty, pars1,
                                      d, 100000)
est_gr2_pars1

est_gr2_pars2 <- estimate_w(loglkl_mvn_penalty, gradient_loglkl_penalty, pars2,
                   d, 100000)
est_gr2_pars2

est_gr2_w <- estimate_w(loglkl_mvn_penalty, gradient_loglkl_penalty, w,
                             d, 100000)
est_gr2_w # initial values as true parameters




# Divide and Conquer ------------------------------------------------------
## get the third pair of initial parameters:
pars3 <- generate_pars_by_range(t = Time, d = d, opt_f = loglkl_mvn_penalty, 
            grad = gradient_loglkl_penalty, maxit = 10000)
pars3

est_gr1_pars3 <- estimate_w(loglkl_mvn_penalty, calc_gradient_num, pars3, d, 100000)
est_gr1_pars3

est_gr2_pars3 <- estimate_w(loglkl_mvn_penalty, gradient_loglkl_penalty, pars3,
                            d, 100000)
est_gr2_pars3




# Tempering Method --------------------------------------------------------
## square the likelihood and take the log, then optimize
est2_gr2_pars1 <- estimate_w(tempering_loglkl_mvn_penalty, 
                             tempering_gradient_loglkl_penalty, 
                             pars1, d, 100000)
est2_gr2_pars1

est2_gr2_pars2 <- estimate_w(tempering_loglkl_mvn_penalty, 
                             tempering_gradient_loglkl_penalty, 
                             pars2, d, 100000)
est2_gr2_pars2

est2_gr2_pars3 <- estimate_w(tempering_loglkl_mvn_penalty, 
                                   tempering_gradient_loglkl_penalty, 
                                   pars3, d, 100000)
est2_gr2_pars3 #they are going crazy

est2_gr2_w <- estimate_w(tempering_loglkl_mvn_penalty, 
                               tempering_gradient_loglkl_penalty, 
                               w, d, 100000)
est2_gr2_w




# Dynamic Linear Models ---------------------------------------------------
est3_gr1_pars1 <- estimate_w(dynamic_loglkl_mvn_penalty, 
                             calc_gradient_num, 
                             pars1, d, 100000)
est3_gr1_pars1




# Parallel programming ----------------------------------------------------
#llply(d, estimate_w, opt_f = loglkl_mvn_penalty, grad = calc_gradient_num, pars = pars3, maxit = 100000, .parallel = 1)




# Assessment --------------------------------------------------------------
mse(est2$par, w)
mse(est4$par, w) 

# difference in estimates:
diff_btw_estimates(est1$par, est2$par)

# plot of differences:
plot_residuals(est1$par, est2$par) 




# MoM (alternative) -------------------------------------------------------
mom_est <- gmm(g = loglkl_mvn_penalty, x = d, t0 = w, 
               gradv = gradient_loglkl_penalty, optfct = "optim", 
               itermax = 1000000)

## did not work as expected




# Plots -------------------------------------------------------------------
## NEW: later
## Concern: there is small discrepancy in estimates, so optim() might be 
## converging to local maxima. To find out, we get 3D plot:

# plot without color
lkl_plot(d) 

# plot w/ color
lkl_plot_col(d)




# Gradients ---------------------------------------------------------------
gr1 <- calc_gradient_num(loglkl_mvn_penalty, w, d, epsilon=10^-6)
gr1

gr2 <- gradient_loglkl_penalty(w, d)
gr2

# little bit discrepancy: is it due the new gradient? which one is causing
# the discrepancy? >> should check it later



















































