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
n <- 40
d_sf_emulator <- data.frame(type = "emulator",
                         total_prcp = rep(seq(from = min(d_sf$total_prcp),
                                              to = max(d_sf$total_prcp),
                                              length.out = n), times = 2)
)

d_sf_emulator <- cbind(d_sf_emulator, 
                       rbind(rep_row(slope_profiles[1,], n),
                             rep_row(slope_profiles[2,], n)))





# MLE estimates -----------------------------------------------------------
# new best
theta_mle_sf <- c(55.06167167, -402.44309188, 
                  -12.67886639,  386.78739175,  0.00197348)

# lkl = 5.490827e+169
#theta_mle <- c(0.011126216, -25.682813936, 4.347803503, 
#               13.044304027, 0.003774865)

# largest lkl = 3.585272e+248
#theta_mle <- c(3.175177e-02, -1.970449e+02,  1.891088e+02,
#              8.627757e+01,  3.819011e-03)

# cvg = 0, lkl = 1.619814e+175
#theta_mle <- c(9.962087e-03, -2.165815e+02, -4.978847e+01,
#               1.983269e+02,  3.487355e-03)

# cvg = 1, lkl = lkl = 2.560958e+248, n = 432
#theta_mle <- c(9.690187e-03, -3.933811e+02, -2.216817e+01,
#               3.863279e+02,  3.750334e-03)

# middle good
#theta_mle <- c(0.010058478, -13.062947983,  -0.684352483,
#               -0.610023781,   0.003783663)

# little good
# theta_mle <- c(0.009203076, -38.852993474,  30.923483013,
#               8.634562892,   0.003552888)

# very bad
#theta_mle <- c(3.929260e-02, -1.940042e+02,  2.535889e+02,
#               2.362047e+00, 3.848922e-03)

# very bad
#theta_mle <- c(0.069210529, -39.963972635, 245.711304322,
#               40.666849998, 0.004035863)

# generate only true parameters (weights + lambda) from theta_mle
basis <- generate_trig_basis(15, 1)
logit_w_sf <- theta_mle_sf[-c(1,5)]%*%basis

# get w
w_sf <- exp(logit_w_sf)/(1+exp(logit_w_sf))

# new theta
theta_sf <- c(theta_mle_sf[1], w_sf, theta_mle_sf[5])





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






# Emulator ----------------------------------------------------------------
## emulate soil loss values
soil_loss <- c()
nhat <- nrow(d_sf_emulator)
for (i in 1:nhat) {
    soil_loss[i] <- emulate_soil_loss(d_sf[,-1], 
                                      d_sf_emulator[i, -1], 
                                      theta_sf,
                                      sigma_sf_inv)
}







# Visualization - did not use ---------------------------------------------
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
dev.copy(pdf,'bwd4_sf.pdf')
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
dev.copy(pdf,'orb4_sf.pdf')
dev.off()







# Add SE ------------------------------------------------------------------
## write a function to find CI borders
ci_sf <- function(d, xhat, mle, sigma_sf_inv, mean_sf) {
    sigma_star <- calc_sigma_star(d, xhat, mle)
    
    # calculate variance
    var <- 1 + mle[17] - sigma_star%*%sigma_sf_inv%*%sigma_star
     
    # get sample of 1000
    ss <- rnorm(1000, mean = mean_sf, sd = sqrt(var))
    
    # get quantiles
    res <- quantile(ss, probs = c(0.025, 0.975), na.rm = TRUE)
    return(res)
}

## add CI for the whole xgrid
nhat <- nrow(d_sf_emulator)
for (i in 1:nhat) {
    ci <- ci_sf(d_sf[,-1], 
                d_sf_emulator[i, -1], 
                theta_sf,
                sigma_sf_inv,
                soil_loss[i])
    d_sf_emulator$q1[i] <- ci[1]
    d_sf_emulator$q2[i] <- ci[2]
}






# Merge datasets into one -------------------------------------------------
d_sf_emulator <- cbind(d_sf_emulator, soil_loss)

# re-standardize precipitaion
d_sf$total_prcp <- exp(d_sf$total_prcp) - 637.56

# re-standardize precipitaion
d_sf_emulator$total_prcp <- exp(d_sf_emulator$total_prcp) - 637.56

# combine all into one
#names(d_sf) <- names(d_sf_emulator) 
#d_sf_mixed <- rbind(d_sf, d_sf_emulator)






# Final Mixed plot --------------------------------------------------------
cols <- c("wepp"="#f04546","emulated"="#17202A","CI"="#62c76b")
# baswood1 hill1 
bswd1h1_mixed <- ggplot() +
    geom_point(data = d_sf[1:12,], 
                aes(x = total_prcp, y = soil_loss, colour = "wepp")) +
    geom_line(data = d_sf_emulator[1:40,],
                aes(x = total_prcp, y = soil_loss, colour = "emulated")) +
    geom_line(data = d_sf_emulator[1:40,],
              aes(x = total_prcp, y = q1, colour = "CI")) +
    geom_line(data = d_sf_emulator[1:40,],
              aes(x = total_prcp, y = q2, colour = "CI")) +
    ggtitle("Soil loss distribution for Basswood1 - hill1 slope profile") + 
    xlab("annual precipitation (mm)") +
    ylab("annual soil loss (kg/m^2)") +
    scale_colour_manual(name="types", values = cols) +
    theme_light()
bswd1h1_mixed

# publish
#dev.copy(pdf,'bwd_sf_mixed.pdf')
dev.copy2pdf(out.type = "pdf")
dev.off()


# orbweaver2 hill3
orb2h3_mixed <- ggplot() +
    geom_point(data = d_sf[13:24,], 
               aes(x = total_prcp, y = soil_loss, colour = "wepp")) +
    geom_line(data = d_sf_emulator[41:80,],
              aes(x = total_prcp, y = soil_loss, colour = "emulated")) +
    geom_line(data = d_sf_emulator[41:80,],
              aes(x = total_prcp, y = q1, colour = "CI")) +
    geom_line(data = d_sf_emulator[41:80,],
              aes(x = total_prcp, y = q2, colour = "CI")) +
    ggtitle("Soil loss distribution for Orbweaver2 - hill3 slope profile") + 
    xlab("annual precipitation (mm)") +
    ylab("annual soil loss (kg/m^2)") +
    scale_colour_manual(name="types", values = cols) +
    theme_light()
orb2h3_mixed

# publish
#dev.copy(pdf,'orb_sf_mixed.pdf')
dev.copy2pdf(out.type = "pdf")
dev.off()

