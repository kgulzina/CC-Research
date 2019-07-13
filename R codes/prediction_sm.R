####### This code plots a GP process with estimated parameters #######
## log is given "README.md" file in github
## Author: stat_cat
## Date: 06/18/19 - present




# Libraries --------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(ggpubr)
library(plotly)
library(lattice)
source("R codes/supplementary_functions.R")
source("R codes/test_on_simulation.R")




# Real data - Scalar case -------------------------------------------------
n <- 50
original_d_sim <- d_simulated
d_sim <- data.frame(type = "simulated",
                   original_d_sim[,-21],
                   y = original_d_sim[,21])

# generate two new new locations: by adding 10^(-6)
d_sim_emulator <- data.frame(type = "emulated",
                             original_d_sim[1:n, -21],
                             id = 1:n)





# MLE estimates -----------------------------------------------------------
# theta = (w_1, w_2, lambda)
theta_mle <- c(-58.0079065,   2.4872092,  32.7025199,   0.0391649)

# generate only true parameters (weights + lambda) from theta_mle
basis <- generate_trig_basis(20, 1)
logit_w_sm <- theta_mle[-4]%*%basis

# get w and v
w_sm <- exp(logit_w_sm)/(1+exp(logit_w_sm))

# new theta
theta_sm <- c(w_sm, theta_mle[4])





# Code predictor - last summer --------------------------------------------
# constants
n <- nrow(d_sim)
sigma_sm <- gcalc_corr_scalar(d_sim[-1], theta_sm[-4]) + 
    diag(x = theta_sm[4], nrow = n, ncol = n)
sigma_sm_inv <- solve(sigma_sm)


# functions
calc_sigma_star <- function(d, xhat, mle) {
    n <- nrow(d)
    y <- ncol(d)
    res <- c()
    for (i in 1:n) {
        res[i] <- exp(-sum(mle[-4]*(d[i,-y] - unlist(xhat))^2))
        # check all vectors if the same
        if (sum(d[i,-y] == unlist(xhat)) == 20) {
            # add nugget
            res[[i]] <- res[[i]] + mle[4]
        }
    }
    return(res)
}



emulate_soil_loss <- function(d, xhat, mle, sigma_sm_inv) {
    # this code emulates soil loss: y*/y >> mean as a prediction
    sigma_star <- calc_sigma_star(d, xhat, mle)
    
    # yhat = sigma_star'sigma^{-1}y
    yhat <- sigma_star%*%sigma_sm_inv%*%d[,21]
    return(yhat)
}


# example
xhat <- d_sim_emulator[1,-c(1,22)]
emulate_soil_loss(d_sim[,-1], xhat, theta_sm, sigma_sm_inv)





# Emulate -----------------------------------------------------------------
## emulator
y <- c()
nhat <- nrow(d_sim_emulator)
for (i in 1:nhat) {
    y[i] <- emulate_soil_loss(d_sim[,-1], 
                                      d_sim_emulator[i, -c(1,22)], 
                                      theta_sm,
                                      sigma_sm_inv)
}

d_sim_emulator <- cbind(d_sim_emulator, y)

# combine all into one
d_sim <- d_sim[1:n, -22] %>% 
    mutate(id = 1:n,
           y = d_sim[1:n, 22])
names(d_sim) <- names(d_sim_emulator) 
d_sim_mixed <- rbind(d_sim, d_sim_emulator)





# Visualization -----------------------------------------------------------
## simulated data: x1, ... , x10
sim <- d_sim %>% 
    ggplot(aes(x = id, y = y)) +
    geom_bar(stat = "identity") +
    xlab("functional input: x values") +
    ylab("simulated y") +
    theme_light()
sim    

# publish
#dev.copy(pdf,'sim.pdf')
#dev.off()


## emulated data: x1, ... , x10
em <- d_sim_emulator %>% 
    ggplot(aes(x = id, y = y)) +
    geom_bar(stat = "identity") +
    xlab("functional input: x values") +
    ylab("emulated y") +
    theme_light()
em

# publish
#dev.copy(pdf,'em.pdf')
#dev.off()

ggarrange(sim, em, nrow = 2)


# scatter plot of mixed as in FF
cols_sm <- c("simulated"="#f04546", "emulated"="#17202A")
sm <- ggplot() +
    geom_point(data = d_sim[1:50,], size = 2.5,
               aes(x = id, y = y, colour = "simulated")) +
    geom_point(data = d_sim_emulator, 
              aes(x = id, y = y, colour = "emulated")) +
    # to remove legend for size
    scale_size(guide = "none") +
    ggtitle("Distribution of the response: y") +
    xlab("(first 50) x values from the simulated data") +
    theme_light()
sm    

# publish
dev.copy2pdf(out.type = "pdf")
dev.off()






