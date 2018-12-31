#### Gradients for optimization ####

# numerical gradient
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

# true gradient for loglkl_with_penalty



































