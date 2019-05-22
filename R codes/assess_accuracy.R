######### This code assesses the estimates for accuracy #########
## Author: stat_cat
## Date: 01/09/19

diff_btw_estimates <- function(est1, est2){
# Gives the mean of squared differences between two estimates. Estimates are
# are different by initial values, gradient, etc.
#
# Args:
#   est1: output of optim() with pars1
#   est2: output of optim() with pars2
#
# Output:
#   msd: mean squared difference (scale)
    
    msd <- mean((est1-est2)^2)
    return(msd)
}




mse <- function(est, w){
# Calculates the MSE of estimates found by particular method. Valid only for
# simulations. 
#
# Args: 
#   est1: output of optim() with pars1
#   w:    simulated weights, in a vector form
#
# Output: 
#   mse: mean squared error (scale)
    
    mse <- mean((est-w)^2)
    return(mse)
    
}




plot_residuals <- function(est, w){
# Plots the residuals = difference between estimates and true weights: w
#
# Args:
#   est1: output of optim() with pars1
#   w:    simulated weights, in a vector form 
#
# Output:
#   graph: residual plot
    
    d <- as.vector(est - w)
    plot(d)
    abline(h = 0, col = "red")
}






