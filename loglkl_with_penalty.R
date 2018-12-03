### Log-likelihood for MVN with penalty ###
library(MASS)
library(mvtnorm)
library(Matrix)

## goals: 
# Optimize functions with loops: remove c() >> Done
# Increase the penalty: rho = 0.99 >> Done
# Set white noise variance: sigmasq = 1 
# Change w(t) from deterministic to sampled
# Remove constants from functions, split variables
# R style guide: make changes
# Split functions into different files
# Derive the real gradient function
# Assess the accuracy of optim with different pars



## conventions: 
# sample size and Time: T << n
n <- 100
Time <- 3 # if 1 -- cannot calculate the sigmasq

# Sigma: covariance matrix for AR(1)
Sigma <- function(t, rho){ # t:vector size; rho is fixed
    sigma <- matrix(NA, nrow = t, ncol = t)
    for(i in 1:t){
        for(j in 1:t){
            sigma[i,j] <- rho^abs(i-j)
        }
    }
    return(sigma)
}

# cov. matrix for MVN simulation: n samples
sigma <- Sigma(Time+1, rho = 0.95)

# simulate
x_values <- mvrnorm(n, rep(0, Time+1), sigma) # vector of length T+1

# Define w(t)=(T-t)/T :straight line with the negative slope: -1/T
# Calculate covariances between simulated X's using correlation form on 
# page3, (1) form.


# a correlation matrix: (using corr(y_x_i, y_x_j)) -- For simulation
calc_corr <- function(val, t, n) { #simulated X's
    omega <- matrix(NA, nrow = n, ncol = n)
    time <- 0:t
    for(i in 1:n) {
        for(j in 1:n) {
            omega[i,j] <- exp(-sum(((t-time+1)/(t+1))*(val[i,]-val[j,])^2))
        }
    }
    return(omega)
}

omega <- calc_corr(x_values, Time, n)

# simulate y
y <- rmvnorm(1, mean = rep(0, n), sigma = omega) 
plot(y[1,]) 
summary(y[1,]) # get the idea



# once we simulate Y's, let's see if they are generated correctly:
# We need positive y_values:
positive_y <- y + abs(min(y))# shift by the minimum value
plot(positive_y[1,])
min1_y <- positive_y + (1-min(positive_y)) # set minimum value of y = 1
plot(min1_y[1,])
# convert to colors:
colors <- colorRampPalette(c("red", "orange", "blue"))(3)[min1_y[1,]]
length(colors) #check if all y have been assigned colors
# add to dataframe:
to_check <- as.data.frame(x_values)
colnames(to_check) <- c("x0","x1")
to_check$colors <- colors

# Plotted are X's (default: n=100) and points are colored by Y:
plot(to_check[,1], to_check[,2], col = to_check$colors)
# close dots should be of the same color



# d: dataset combined: functional input and scale output
d <- cbind(x_values,t(y)) #colnames(d) <- c("t0", "t1", "y")

# w(t)=rho*w(t-1) + z_t, z_t ~ N(0,sigma_square)
# to find sigma_square:
find_sigmasq <- function(rho, t){
    time <- seq(0, t, by = 1)
    v <- (t-time+1)/(t+1)
    n <- length(v)
    z <- numeric(n-1) 
    for(i in 1:(n-1)) {
        z[i] <- v[i+1]-rho*v[i]
    }
    return(var(z))
}

find_sigmasq(0.99, 2)

# need omega -- general w/ w(t) as a vector. 
gcalc_corr <- function(d,w) { #retriev x's from d
    n <- nrow(d)
    time <- ncol(d)-1
    omega <- matrix(NA, nrow = n, ncol = n)
    
    # remove last column: y
    for(i in 1:n){
        for(j in 1:n){
            omega[i,j] <- exp(-sum(w*(d[i,-ncol(d)]-d[j,-ncol(d)])^2))
        }
    }
    return(omega)
}

pars <- seq(1, 1/(Time+1), len = Time+1)
trial <- gcalc_corr(d, pars) 
rankMatrix(trial)



# new log_lkl w/ penalty term to optimize:
loglkl_mvn_penalty <- function(w,d) { # 
    
    time <- ncol(d)-2
    # calculate the covariance matrix for yw
    omega <- gcalc_corr(d,w)
    # find sigma_square for the prior covariance matrix
    ssq <- find_sigmasq(0.99, time)
    # calculate the covariance matrix for w
    sigma <- ssq/(0.36)*Sigma(time+1, 0.99) ## remove constants
    
    # data model: log_likelihood 
    p_yw <- dmvnorm(d[,ncol(d)], mean = rep(0, nrow(d)), sigma = omega, log = TRUE)  
    
    # prior on w: log_likelihood
    p_w <- dmvnorm(w, mean = rep(0, time+1), sigma = sigma , log = TRUE)
    
    # the posterior which will be maximized
    result <- p_yw + p_w 
    
    return(result)
}


# w = (w_0, w_1, ... , w_100)
# set initial values for the parameters
pars <- seq(1, 1/(Time+1), len = Time+1)
pars2 <- seq(1,0.001, len = Time+1) #gives approximately the same results!

# use log_likelihood w/ penalty:
opt <- optim(par = pars, loglkl_mvn_penalty, d = d,
              control = list(fnscale = -1, maxit=4000))
opt$convergence 
opt$par
pars


# compare w/ pars2:
opt2 <- optim(par = pars2, loglkl_mvn_penalty, d = d,
             control = list(fnscale = -1, maxit=4000))
opt2$convergence 
opt2$par
pars2



# add numerical gradient to optim() 
calc_gradient_num <- function(w,d,epsilon=10^-8){
    # calculates the gradient numerically
    n <- length(w)
    gr <- numeric(n) 
    for(i in 1:n) {
        h <- rep(0,n); h[i] <- epsilon
        gr[i] <- (loglkl_mvn_penalty(w+h,d)-
                      loglkl_mvn_penalty(w,d))/epsilon
    }
    return(gr)
}


######## Automation ##########

# constant terms: simuated data
simulate_d <- function(Time, n){
    
    sigma <- Sigma(Time+1, rho = 0.95)
    
    # simulate X's
    x_values <- mvrnorm(n, rep(0, Time+1), sigma)
    
    #calculate omega for Y
    omega <- calc_corr(x_values, Time, n)
    
    # simulate
    y <- rmvnorm(1, mean = rep(0, n), sigma = omega)
    
    # combine datasets: functional input and scale output
    d <- cbind(x_values,t(y))
    
    return(d)
    
}

d <- simulate_d(10, 100)


# constant term: pars
pars <- seq(1, 1/(Time+1), len = Time+1) # true ones


estimate_w <- function(opt_f, grad, pars, d) { 
    
    # returns optim() output with tru parameter values
    opt <- optim(par = pars, opt_f, d = d, control = list(fnscale = -1,
                                                     maxit=10000),
                 gr = grad) # added gradient function
    
    print(pars)
    return(opt)
}

estimate_w(loglkl_mvn_penalty, calc_gradient)

































