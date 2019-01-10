


library(plotly)
library(lattice)


d <- cbind(x_values,t(y))
#colnames(d) <- c("t0", "t1", "y")

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
                w[t+1] <- 0.8*w[t] + rnorm(1, w[t], .001)
                ss <- ss + w[t]*(d[i,t]-d[j,t])^2
            }
            sigma[i,j] <- exp(-ss)
        }
    }
    
    result <- -(1/2)*(n*log(2*pi)+determinant(sigma, logarithm = TRUE)$modulus) + 
        (t(d[,t+1])%*%ginv(sigma)%*%d[,t+1])/2
    # used ginv()
}




# to see the estimates using any likelihood function for range of Times:
w_estimate <- function(Time, n, opt_f) {
    sigma <- Sigma(Time+1, rho = 0.95)
    
    # simulate
    x_values <- rmvnorm(n, mean = rep(0, Time+1), sigma = sigma)
    
    #calculate omega for Y
    omega <- calc_corr(x_values, Time, n)
    
    # simulate
    y <- rmvnorm(1, mean = rep(0, n), sigma = omega)
    
    # combine datasets: functional input and scale output
    d <- cbind(x_values,t(y))
    
    # set initial values for the parameters
    pars <- seq(1, 1/(Time+1), len = Time+1)
    
    # find gradient
    #gr <- calc_gradient(opt_f, 10^-10, d)
    
    # optimize
    opt <- optim(par = pars, opt_f, d = d, control = list(fnscale = -1,
                maxit=1000, trace=1), method = "L-BFGS-B",
                gr = calc_gradient, lower = rep(-1, Time+1), 
                upper = rep(2, Time+1))
    
    print(pars)
    return(opt)
}

w_estimate(3, 100, loglkl_mvn_penalty)



### sketch the graph of a lkl function with penalty term: 
# general p-dimensional code 
# keeps selected axes constant 


# find values of lkl function:
lkl_values <- function(d){
    seq <- seq(0, 1, 0.01)
    z <- matrix(NA, length(seq), length(seq))
    for(i in 1:length(seq)) {
        for(j in 1:length(seq)) {
            z[j,i] <- loglkl_mvn_penalty(c(seq[i], seq[j]), d)
        }
    }
    return(z) # z should be matrix for persp()
}

lkl_values(d)


# special lkl function for plotting:
loglkl_mvn_penalty <- function(w1,w2,d) {
    n <- nrow(d)
    time <- ncol(d)-1
    sigma <- matrix(NA, ncol = n, nrow = n)
    s <- find_sigma(0.8, time-1)
    
    for(i in 1:n){
        for(j in 1:n){
                ss <- w1*(d[i,1]-d[j,1])^2 + w2*(d[i,2]-d[j,2])^2
            }
            sigma[i,j] <- exp(-ss)
        }
    
    result <- -(1/2)*(n*log(2*pi)+determinant(sigma, logarithm = TRUE)$modulus) + 
        (t(d[,t+1])%*%ginv(sigma)%*%d[,t+1])/2
    # used ginv()
}


# plot w/out color
lkl_plot <- function(Time,n) {
    sigma <- Sigma(Time+1, rho = 0.95)
    
    # simulate
    x_values <- mvrnorm(n, rep(0, Time+1), sigma)
    
    #calculate omega for Y
    omega <- calc_corr(x_values, Time, n)
    
    # simulate
    y <- rmvnorm(1, mean = rep(0, n), sigma = omega)
    
    # combine datasets: functional input and scale output
    d <- cbind(x_values,t(y))
    
    # keep one axis constant, calculate z from outside
    z <- lkl_values(d)
    
    persp(seq(0, 1, 0.01), seq(0, 1, 0.01), z, phi = 45, theta = 0,
          xlab = "weight at time i", ylab = "weight at time j",
          main = "Log-likelihood with penalty", zlim = range(z, na.rm = 0)
          # plot -z, since lkl is negative
    )
}

lkl_plot(1,50) # non-differentiable surface (very bad)!!!!!!


# plot w/ color
lkl_plot_col <- function(Time,n) {
    sigma <- Sigma(Time+1, rho = 0.95)
    
    # simulate
    x_values <- mvrnorm(n, rep(0, Time+1), sigma)
    
    #calculate omega for Y
    omega <- calc_corr(x_values, Time, n)
    
    # simulate
    y <- rmvnorm(1, mean = rep(0, n), sigma = omega)
    
    # combine datasets: functional input and scale output
    d <- cbind(x_values,t(y))
    
    # keep one axis constant, calculate z from outside
    z <- lkl_values(d)
    
    wireframe(seq(0, 1, 0.01), seq(0, 1, 0.01), z, phi = 45, theta = 45,
          xlab = "weight at time i", ylab = "weight at time j",
          main = "Log-likelihood with penalty", colorkey = TRUE)
}

lkl_plot_col(1,50)















































