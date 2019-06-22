####### This code plots a GP process with estimated parameters #######
## log is given "README.md" file in github
## Author: stat_cat
## Date: 06/18/19 - present




# Libraries --------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(plotly)
library(lattice)
source("R codes/supplementary_functions.R")




# Real data - Scalar case -------------------------------------------------
d <- read.csv("data/final_csv/annual_soil_loss_not_stdrzd.csv")
d_sf <- data.frame(type = "wepp",
                   total_prcp = apply(d[,5:369], 1, sum),
                   d[,370:384],
                   soil_loss = d[,4])

# standardize precipitation amount
d_sf[,2] <- log(d_sf[,2] + min(d_sf[d_sf[,2] != 0,2]))

# choose two differenet slope profiles: Basswood1&hill1 and Orbweaver2&hill3
slope_profiles <- unique(d[,370:384])[c(1,33),]

# add new slope and precipitation grid
d_sf_emulator <- data.frame(type = "emulator",
                            total_prcp = rep(seq(from = min(d_sf$total_prcp),
                                                 to = max(d_sf$total_prcp),
                                                 length.out = 12), times = 2)
)

d_sf_emulator <- cbind(d_sf_emulator, 
                       rbind(rep_row(slope_profiles[1,], 12),
                             rep_row(slope_profiles[2,], 12)))





# MLE estimates -----------------------------------------------------------
# theta = (w_1, w_2, lambda)
theta_mle <- c(0.011126216, -25.682813936, 4.347803503, 
               13.044304027, 0.003774865)

# generate only true parameters (weights + lambda) from theta_mle
basis <- generate_trig_basis(15, 1)
logit_w_sf <- theta_mle[-c(1,5)]%*%basis

# get w
w_sf <- exp(logit_w_sf)/(1+exp(logit_w_sf))

# new theta
theta_sf <- c(theta_mle[1], w_sf, theta_mle[5])





# Code predictor - last summer --------------------------------------------
# constants
n <- nrow(d_sf)
sigma_sf <- gcalc_corr_scalar(d_sf[,-1], theta_sf[-17]) + 
    diag(x = theta_sf[17], nrow = n, ncol = n)
sigma_sf_inv <- solve(sigma_sf)


# functions
calc_sigma_star <- function(d, xhat, mle) {
    n <- nrow(d)
    y <- ncol(d)
    res <- c()
    for (i in 1:n) {
        res[i] <- exp(-sum(mle[-17]*(d[i,-y] - unlist(xhat))^2))
        # check all vectors if the same
        if (sum(d[i,-y] == unlist(xhat)) == 16) {
            # add nugget
            res[[i]] <- res[[i]] + mle[17]
        }
    }
    return(res)
}



emulate_soil_loss <- function(d, xhat, mle, sigma_sf_inv) {
    # this code emulates soil loss: y*/y >> mean as a prediction
    sigma_star <- calc_sigma_star(d, xhat, mle)
    
    # yhat = sigma_star'sigma^{-1}y
    yhat <- sigma_star%*%sigma_sf_inv%*%d[,17]
    return(yhat)
}


# example
xhat <- d_sf_emulator[1,-1]
emulate_soil_loss(d_sf[,-1], xhat, theta_sf, sigma_sf_inv)





# Emulate -----------------------------------------------------------------
## emulator
soil_loss <- c()
nhat <- nrow(d_sf_emulator)
for (i in 1:nhat) {
    soil_loss[i] <- emulate_soil_loss(d_sf[,-1], 
                                      d_sf_emulator[i, -1], 
                                      theta_sf,
                                      sigma_sf_inv)
}

d_sf_emulator <- cbind(d_sf_emulator, soil_loss)

# restandardize precipitaion
d_sf$total_prcp <- exp(d_sf$total_prcp) - 637.56

# standardize precipitaion
d_sf_emulator$total_prcp <- exp(d_sf_emulator$total_prcp) - 637.56

# combine all into one
names(d_sf) <- names(d_sf_emulator) 
d_sf_mixed <- rbind(d_sf, d_sf_emulator)








# Visualization -----------------------------------------------------------
## wepp data: baswood1 & hill1
wepp_bsw1 <- d_sf[1:12,-1] %>% ggplot(aes(x = total_prcp, y = soil_loss)) +
    geom_smooth() + #alternatively use geom_line() or geom_point() 
    ggtitle("Soil loss distribution for Basswood1 - hill1 slope profile") + 
    xlab("annual precipitation (mm)") +
    ylab("annual soil loss (kg/m^2)")
wepp_bsw1

## wepp data: orbweaver2 & hill3
wepp_orb3 <- d_sf[385:396,-1] %>% ggplot(aes(x = total_prcp, y = soil_loss)) +
    geom_smooth() +
    ggtitle("Soil loss distribution for Orbweaver2 - hill3 slope profile") + 
    xlab("annual precipitation (mm)") +
    ylab("annual soil loss (kg/m^2)")
wepp_orb3

## emulator data: baswood1 & hill1
em_bsw1 <- d_sf_emulator[1:12,-1] %>% ggplot(aes(x = total_prcp, y = soil_loss)) +
    geom_smooth() + 
    ggtitle("Soil loss distribution for Basswood1 - hill1 slope profile") + 
    xlab("annual precipitation (mm)") +
    ylab("annual soil loss (kg/m^2)")
em_bsw1

## emulator data: orbweaver2 & hill3
em_orb3 <- d_sf_emulator[13:24,-1] %>% ggplot(aes(x = total_prcp, y = soil_loss)) +
    geom_smooth() +
    ggtitle("Soil loss distribution for Orbweaver2 - hill3 slope profile") + 
    xlab("annual precipitation (mm)") +
    ylab("annual soil loss (kg/m^2)")
em_orb3


## plot mixed: baswood1 & hill1
d_sf_mixed[c(1:12, 433:444),] %>% ggplot(aes(x = total_prcp, 
                                             y = soil_loss,
                                             color = type)) +
    geom_smooth() + #alternatively use geom_line() or geom_point() 
    ggtitle("Soil loss distribution for Basswood1 - hill1 slope profile") + 
    xlab("annual precipitation (mm)") +
    ylab("annual soil loss (kg/m^2)") +
    theme_linedraw()

# publish
dev.copy(pdf,'bwd1_sf.pdf')
dev.off()


## plot mixed: orbweaver2 & hill3
d_sf_mixed[c(385:396, 445:456),] %>% ggplot(aes(x = total_prcp, 
                                                y = soil_loss,
                                                color = type)) +
    geom_smooth() + #alternatively use geom_line() or geom_point() 
    ggtitle("Soil loss distribution for Orbweaver2 - hill3 slope profile") + 
    xlab("annual precipitation (mm)") +
    ylab("annual soil loss (kg/m^2)") +
    theme_linedraw()

# publish
dev.copy(pdf,'prcp.pdf')
dev.off()











