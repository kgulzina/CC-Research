# PART1
# Simulate X = {x_1, x_2, ... , x_n} where each x_i is a vector of length T+1,
# s.t length(x_i)=T+1. X ~ GP(0, Sigma). Sigma is AR(1) with rho=0.95

library(MASS) #github
library(mvtnorm)
library(Matrix)

# a covariance matrix

Sigma <- function(t, rho){ # n:sample size; rho is fixed
    sigma <- matrix(NA, nrow = t, ncol = t)
    for(i in 1:t){
        for(j in 1:t){
            sigma[i,j] <- rho^abs(i-j)
        }
    }
    return(sigma)
}


# sample size and Time: T << n
n <- 100
Time <- 3
    
# sigma for MVN simulation: n samples
sigma <- Sigma(Time+1, rho = 0.95)

# simulate
x_values <- mvrnorm(n, rep(0, Time+1), sigma) # vector of length
## T+1. I dont know if its correct since the sample covariances are negative 
## for some. Should be close to 0: 0.95^n-1 --->0 as n-->inf

# try another function for simulation:
x_values2 <- rmvnorm(n, mean = rep(0, Time+1), sigma = sigma)
## it seems that there is no difference between two simulations





#PART2
# Define w(t)=(T-t)/T :straight line with the negative slope: -1/T
# Calculate covariances between simulated X's using correlation form on 
# page3, (1) form. 

# a covarince matrix using d(x_i,x_j): distance only

calc_dist <- function(val, t, n) { #simulated X's
    omega <- matrix(NA, nrow = n, ncol = n)
    time <- 0:t
    for(i in 1:n) {
        for(j in 1:n) {
            if(i!=j)
            omega[i,j] <- sum(((t-time)/t)*(val[i,]-val[j,])^2)
            else omega[i,j] <- 1 # confusion???
        }
    }
    return(omega)
}

dist_mtx <- calc_dist(x_values, Time, n)   



# a correlation matrix: (using corr(y_x_i, y_x_j))
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

rankMatrix(omega) # is FULL rank only when T>1!! (WHY?)


#PART3
y <- rmvnorm(1, mean = rep(0, n), sigma = omega) 
plot(y[1,])
summary(y[1,])

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

## trying the way 2:
abs_y <- abs(y) + (1-min(abs(y)))
plot(abs_y[1,])
colors2 <- colorRampPalette(c("red", "orange", "blue"))(3)[abs_y[1,]]
length(colors2) #check if all y have been assigned colors
# add to dataframe:
to_check2 <- as.data.frame(x_values)
colnames(to_check2) <- c("x0","x1")
to_check2$colors <- colors2
plot(to_check[,1], to_check[,2], col = to_check$colors)
# close dots appear to be of the same color depending on a cluster they
# form within the plane. 
### Correct? Dr. Niemi 

#PART4
# Find MLE of w(t), t = 0:100
# Write definition of Y~GP(0,Omega) in general as in summmer. Use optim.

# combine datasets: functional input and scale output
d <- cbind(x_values,t(y))
#colnames(d) <- c("t0", "t1", "y")

# Likelihood function for the multivariate normal: Y
lkl_mvn <- function(w,d) {
    n <- nrow(d)
    t <- ncol(d)-1
    sigma <- matrix(NA, ncol = n, nrow = n)
    for(i in 1:n){
        for(j in 1:n){
            ss <- 0
            for(t in 1:t){
                ss <- ss + w[t]*(d[i,t]-d[j,t])^2
            }
            sigma[i,j] <- exp(-ss)
        }
    }
    
    result <- sqrt(1/((2*pi)^n)*det(sigma))*exp((t(d[,t+1])%*%ginv(sigma)%*%d[,t+1])/2) 
    # solve(sigma) gives error, used ginv() ??? is it valid?
    return(result)
}

# Log_likelihood function for the multivariate normal: Y
loglkl_mvn <- function(w,d) {
    n <- nrow(d)
    t <- ncol(d)-1
    sigma <- matrix(NA, ncol = n, nrow = n)
    for(i in 1:n){
        for(j in 1:n){
            ss <- 0
            for(t in 1:t){
                ss <- ss + w[t]*(d[i,t]-d[j,t])^2
            }
            sigma[i,j] <- exp(-ss)
        }
    }
    # write the restriction on w(t) !!!
    result <- -(1/2)*(n*log(2*pi)+determinant(sigma, logarithm = TRUE)$modulus) + 
        (t(d[,t+1])%*%ginv(sigma)%*%d[,t+1])/2
    # used ginv()
}


# w = (w_0, w_1, ... , w_100)
# set initial values for the parameters
pars <- seq(1, 1/(Time+1), len = Time+1)
    
# use optimize to Maximize
opt <- optim(par = pars, lkl_mvn, d = d, control = list(fnscale = -1))
# since optim() minimizes as a default, we will set control fnscale to -1. 
mle <- opt$par
mle

# use log_likelihood:
opt2 <- optim(par = pars, loglkl_mvn, d = d, control = list(fnscale = -1))

mle2 <- opt2$par
opt2$convergence #didn't converge
mle2  #different and close to the true values


#10/31/18
### As we see optimization did not converge correctly, so we will try to 
### optimize new loss function which has some penalty term. It should do 
### regularization on the estimation of the parameters. 

# As a penalty we will assign a prior distribution to the weights w(t): 
# Autoregressive model: AR(1), with coefficient 0<rho<1 s.t rho=0.8. We have
# to find sigma_hat_square based on the true w(t) and rho=0.8

# w(t)=rho*w(t-1) + delta_t, delta_t~N(0,sigma_square)
# to find sigma_square:
find_sigma <- function(rho, t){
    time <- seq(0, t, by = 1)
    v <- (t-time+1)/(t+1)
    for(i in 1:(length(v)-1)) {
        delta <- v[i+1]-rho*v[i]
    }
    return(sd(v))
}

find_sigma(0.8, Time)

# new log_lkl w/ penalty term to optimize:
loglkl_mvn_penalty <- function(w,d) {
    n <- nrow(d)
    time <- ncol(d)-1
    sigma <- matrix(NA, ncol = n, nrow = n)
    s <- find_sigma(0.8, time-1)
    
    for(i in 1:n){
        for(j in 1:n){
            ss <- 0
            for(t in 1:time){
                #restriction on w(t)
                w[t+1] <- 0.8*w[t] + rnorm(1, 0, s)
                ss <- ss + w[t]*(d[i,t]-d[j,t])^2
            }
            sigma[i,j] <- exp(-ss)
        }
    }
    
    result <- -(1/2)*(n*log(2*pi)+determinant(sigma, logarithm = TRUE)$modulus) + 
        (t(d[,t+1])%*%ginv(sigma)%*%d[,t+1])/2
    # used ginv()
}

# use log_likelihood w/ penalty:
opt3 <- optim(par = pars, loglkl_mvn_penalty, d = d, control = list(fnscale = -1))
opt3$convergence # didn't converge
opt3$par # pretty like close 
pars


# to see the estimates using any likelihood function for range of Times:
w_estimate <- function(Time, n, opt_f) {
    sigma <- Sigma(Time+1, rho = 0.95)
    
    # simulate
    x_values <- mvrnorm(n, rep(0, Time+1), sigma)
    
    #calculate omega for Y
    omega <- calc_corr(x_values, Time, n)
    
    # simulate
    y <- rmvnorm(1, mean = rep(0, n), sigma = omega)
    
    # combine datasets: functional input and scale output
    d <- cbind(x_values,t(y))
    
    # set initial values for the parameters
    pars <- seq(1, 1/(Time+1), len = Time+1)
    
    # optimize
    opt <- optim(par = pars, opt_f, d = d, control = list(fnscale = -1,
                                                          maxit=1000,
                                                          trace=1))
    
    print(opt$convergence)
    print(pars)
    return(opt)
}

w_estimate(10, 100, loglkl_mvn_penalty)


### sketch the graph of a lkl function with penalty term: 
# general p-dimensional code 
# keeps selected axes constant 















































