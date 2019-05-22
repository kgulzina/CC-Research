######### This code plot 3D graph of log likelihood function #########
## Author: stat_cat
## Date: 01/09/19

## This function works Time = 1 (only)!


lkl_values <- function(d){
# Gives the set of output values of log-likelihood function
#
# Args:
#   d: 
# Output:
#   z: matrix of function values
    
    seq <- seq(0, 1, 0.001)
    z <- matrix(NA, length(seq), length(seq))
    for(i in 1:length(seq)) {
        for(j in 1:length(seq)) {
            z[j,i] <- loglkl_mvn_penalty(c(seq[i], seq[j]), d)
        }
    }
    return(z) 
}




lkl_plot <- function(d) {
# Plots black-white 3D graph of log-likelihood function
#
# Args:
#   d: data frame, last column is response Y, others are input X's
# Output:
#   graph: 3D plot 
    
    # keep one axis constant, calculate z from outside
    z <- lkl_values(d)
    
    persp(seq(0, 1, 0.001), seq(0, 1, 0.001), z, phi = 45, theta = 45,
          xlab = "weight at time i", ylab = "weight at time j",
          main = "Log-likelihood with penalty", zlim = c(0, 10e+10))
}




lkl_plot_col <- function(Time,n) {
# Plots colorful 3D graph of log-likelihood function
#
# Args:
#   d: data frame, last column is response Y, others are input X's
# Output:
#   graph: 3D plot 
    
    # keep one axis constant, calculate z from outside
    z <- lkl_values(d)
    
    wireframe(seq(0, 1, 0.001), seq(0, 1, 0.001), z, phi = 45, theta = 45,
          xlab = "weight at time i", ylab = "weight at time j",
          main = "Log-likelihood with penalty", zlim = c(0, 10e+10),
          colorkey = TRUE)
}
















































