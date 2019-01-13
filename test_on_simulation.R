####### This code tests "loglkl_with_penalty" with the simulated data #######
## log is given "README.md" file in github
## Author: stat_cat
## Date: 09/03/18 - present


library(MASS)
library(mvtnorm)
library(Matrix)
library(tmvtnorm)
library(plotly)
library(lattice)


## conventions: 
## sample size and Time: T << n
n <- 100
Time <- 5 
rho <- 0.99
sigmasq <- 1

## to use supplementary functions:
source('supplementary_functions.R')
source('loglkl_with_penalty.R')
source("assess_accuracy.R")
source("3D_plot_loglkl.R")


## simulate weights:
w <- simulate_w_trnctd(Time)

## Approach 1 to simulate w:
w1 <- simulate_log_w_stan(Time)

## Approach 2 to simulate w:
w <- simulate_log_w(Time) #current

## simulate data:
d <- simulate_d(Time, n, w)




#### Some manipulations to see if Y's are sampled correctly ####
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




## set initial values:
pars <- seq(1, 1/(Time+1), len = Time+1)
pars

## compare two estimates with different initial values:
est1 <- estimate_w(loglkl_mvn_penalty, calc_gradient_num, pars, d, 100000)
est1
est2 <- estimate_w(loglkl_mvn_penalty, calc_gradient_num, w, d, 100000)
est2



## compare two estimates with different gradient and pars:
est3 <- estimate_w(loglkl_mvn_penalty, gradient_loglkl_penalty, pars,
                                      d, 100000)
est3
est4 <- estimate_w(loglkl_mvn_penalty, gradient_loglkl_penalty, w,
                   d, 100000)
est4


### mse for different gradients
mse(est2$par, w)
mse(est4$par, w) ### No difference??? are you kidding? !!!! There is 
## improvement in the runtime!!!!




# difference in estimates:
diff_btw_estimates(est1$par, est2$par)

# plot of differences:
plot_residuals(est1$par, est2$par) 





## NEW: later
## Concern: there is small discrepancy in estimates, so optim() might be 
## converging to local maxima. To find out, we get 3D plot:

# plot without color
lkl_plot(d) 

# plot w/ color
lkl_plot_col(d)


### check the gradients:
gr1 <- calc_gradient_num(loglkl_mvn_penalty, w, d, epsilon=10^-6)
gr1

gr2 <- gradient_loglkl_penalty(w, d)
gr2





























































