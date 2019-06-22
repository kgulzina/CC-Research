######### This code contains WEPP data related information #########
## Author: stat_cat
## Date: 06/07/19





# Libraries ---------------------------------------------------------------
library(lubridate)
library(ggplot2)
library(ggthemes)
library(ggpubr)



# Read data ---------------------------------------------------------------
source("R codes/wepp_data.R")



# Visualize covariates ----------------------------------------------------
## add new column = days' order
climate_for_visual <- climate[-leap,] %>% 
    mutate(days = rep(1:365, 12))

## standardize prcp if needed: Dr. Niemi ????
# standardize prcp: by taking log(prcp + smallest non-zero prcp)
for (i in 2007:2018) {
    rows <- which(climate_for_visual$year == i)
    climate_for_visual$prcp[rows] <- 
        log(climate_for_visual$prcp[rows] + 
                min(climate_for_visual$prcp[rows][climate_for_visual$prcp[rows] != 0]))
}

## plot two years: 2007 and 2018
prcp_2007 <- climate_for_visual %>% filter(year == 2007) %>% 
    select(days, prcp) %>% 
    ggplot(aes(x = days, y = prcp)) + 
    geom_line() + ggtitle("Daily precipitation in 2007") +
    xlab("consecutive days") + ylab("precipitation (mm)") + 
    theme_light()

prcp_2018 <- climate_for_visual %>% filter(year == 2018) %>% 
    select(days, prcp) %>% 
    ggplot(aes(x = days, y = prcp)) + 
    geom_line() + ggtitle("Daily precipitation in 2018") +
    xlab("consecutive days") + ylab("precipitation (mm)") + 
    theme_linedraw()

# combine into one
ggarrange(prcp_2007, prcp_2018, nrow = 2)

# publish
#dev.copy(pdf,'prcp.pdf')
#dev.off()


## plot two hillslope slope profiles: Basswood1 hill1, Orbweaver2 hill3
slope1 <- num_slope[1, -c(1, 2)] %>% t() %>% 
    as.data.frame() %>% 
    mutate(n = 1:15)
slope2 <- num_slope[33, -c(1, 2)] %>% t() %>% 
    as.data.frame() %>% 
    mutate(n = 1:15)

# set correct column names
colnames(slope1) <- c("slope", "n")
colnames(slope2) <- c("slope", "n")

# ggplots
bswd1h1 <- slope1 %>% ggplot(aes(x = n, y = slope)) + 
    geom_line() + ggtitle("Slope profile of Basswood1 hill1") +
    xlab("n") + ylab("slope steepness") + 
    theme_light()

orbw2h3 <- slope2 %>% ggplot(aes(x = n, y = slope)) + 
    geom_line() + ggtitle("Slope profile of Orbweaver2 hill3") +
    xlab("n") + ylab("slope steepness") + 
    theme_light()

# combine into one
ggarrange(bswd1h1, orbw2h3, nrow = 2)

# publish
dev.copy(pdf,'slope.pdf')
dev.off()




# Visualize output --------------------------------------------------------
## there are 432 soil loss values, for 36 hillslopes and 12 years.
## plot two hillslopes: Basswood1 hill1 and Orbweaver2 hill3
# ggplots
soil1 <- soil_loss %>% filter(watershed == "Basswood1",
                              hill == "hill1") %>% 
    select(soil_loss, year) %>% 
    ggplot(aes(x = year, y = soil_loss)) + 
    geom_point(size = 2) + ggtitle("Soil loss for Basswood1 hill1") +
    xlab("year") + ylab("soil loss (kg/m^2)") + 
    theme_light()

soil2 <- soil_loss %>% filter(watershed == "Orbweaver2",
                              hill == "hill3") %>% 
    select(soil_loss, year) %>% 
    ggplot(aes(x = year, y = soil_loss)) + 
    geom_point(size = 2) + ggtitle("Soil loss for Orbweaver2 hill3") +
    xlab("year") + ylab("soil loss (kg/m^2)") + 
    theme_light()

# combine into one
ggarrange(soil1, soil2, nrow = 2)

# publish
dev.copy(pdf,'soil_loss.pdf')
dev.off()


## soil_loss all in all
soil_loss %>% ggplot(aes(y = soil_loss, x = year)) + 
    geom_boxplot() + ggtitle("Annual soil loss amount for 36 hillslopes") +
    xlab("year") + ylab("soil loss (kg/m^2)") + 
    theme_light()

# publish
dev.copy(pdf,'soil_loss_all.pdf')
dev.off()





# Visualize hillslope by height and slope ---------------------------------































