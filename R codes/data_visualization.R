######### This code contains WEPP data related information #########
## Author: stat_cat
## Date: 06/07/19





# Libraries ---------------------------------------------------------------
library(lubridate)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(dplyr)



# Read data ---------------------------------------------------------------
source("R codes/wepp_data.R")



# Visualize precipitation -------------------------------------------------
## add new column = days' order
climate_for_visual <- climate[-leap,] %>% 
    mutate(days = rep(1:365, 12))

## standardize prcp if needed: Dr. Niemi (don't need!)
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
    geom_point(col = "blue") + 
    geom_line(aes(x = days, y = prcp, alpha = 0.001), col = "#73DCFF") +
    ggtitle("Daily precipitation in 2007") +
    xlab("") + ylab("") + 
    theme_light() + scale_alpha(guide = "none")
prcp_2007

prcp_2018 <- climate_for_visual %>% filter(year == 2018) %>% 
    select(days, prcp) %>% 
    ggplot(aes(x = days, y = prcp)) + 
    geom_point(col = "blue") + 
    geom_line(aes(x = days, y = prcp, alpha = 0.001), col = "#73DCFF") +
    ggtitle("Daily precipitation in 2018") +
    xlab("") + ylab("") +
    theme_light() + scale_alpha(guide = "none")

# combine into one
ggarrange(prcp_2007, prcp_2018, nrow = 2)

# publish
#dev.copy(pdf,'prcp.pdf')
dev.copy2pdf(out.type = "pdf")
dev.off()

### ADD ALL YEARS : Write a function here later!!!!!
prcp_2008 <- climate_for_visual %>% filter(year == 2008) %>% 
    select(days, prcp) %>% 
    ggplot(aes(x = days, y = prcp)) + 
    geom_point(col = "blue") + 
    geom_line(aes(x = days, y = prcp, alpha = 0.001), col = "#73DCFF") +
    ggtitle("Daily precipitation in 2008") +
    xlab("") + ylab("") +
    theme_light() + scale_alpha(guide = "none")

prcp_2009 <- climate_for_visual %>% filter(year == 2009) %>% 
    select(days, prcp) %>% 
    ggplot(aes(x = days, y = prcp)) + 
    geom_point(col = "blue") + 
    geom_line(aes(x = days, y = prcp, alpha = 0.001), col = "#73DCFF") +
    ggtitle("Daily precipitation in 2009") +
    xlab("") + ylab("") + 
    theme_light() + scale_alpha(guide = "none")

prcp_2010 <- climate_for_visual %>% filter(year == 2010) %>% 
    select(days, prcp) %>% 
    ggplot(aes(x = days, y = prcp)) + 
    geom_point(col = "blue") + 
    geom_line(aes(x = days, y = prcp, alpha = 0.001), col = "#73DCFF") +
    ggtitle("Daily precipitation in 2010") +
    xlab("") + ylab("") + 
    theme_light() + scale_alpha(guide = "none")

prcp_2011 <- climate_for_visual %>% filter(year == 2011) %>% 
    select(days, prcp) %>% 
    ggplot(aes(x = days, y = prcp)) + 
    geom_point(col = "blue") + 
    geom_line(aes(x = days, y = prcp, alpha = 0.001), col = "#73DCFF") +
    ggtitle("Daily precipitation in 2011") +
    xlab("") + ylab("") + 
    theme_light() + scale_alpha(guide = "none")

prcp_2012 <- climate_for_visual %>% filter(year == 2012) %>% 
    select(days, prcp) %>% 
    ggplot(aes(x = days, y = prcp)) + 
    geom_point(col = "blue") + 
    geom_line(aes(x = days, y = prcp, alpha = 0.001), col = "#73DCFF") +
    ggtitle("Daily precipitation in 2012") +
    xlab("") + ylab("") + 
    theme_light() + scale_alpha(guide = "none")

prcp_2013 <- climate_for_visual %>% filter(year == 2013) %>% 
    select(days, prcp) %>% 
    ggplot(aes(x = days, y = prcp)) + 
    geom_point(col = "blue") + 
    geom_line(aes(x = days, y = prcp, alpha = 0.001), col = "#73DCFF") +
    ggtitle("Daily precipitation in 2013") +
    xlab("") + ylab("") + 
    theme_light() + scale_alpha(guide = "none")

prcp_2014 <- climate_for_visual %>% filter(year == 2014) %>% 
    select(days, prcp) %>% 
    ggplot(aes(x = days, y = prcp)) + 
    geom_point(col = "blue") + 
    geom_line(aes(x = days, y = prcp, alpha = 0.001), col = "#73DCFF") +
    ggtitle("Daily precipitation in 2014") +
    xlab("") + ylab("") + 
    theme_light() + scale_alpha(guide = "none")

prcp_2015 <- climate_for_visual %>% filter(year == 2015) %>% 
    select(days, prcp) %>% 
    ggplot(aes(x = days, y = prcp)) + 
    geom_point(col = "blue") + 
    geom_line(aes(x = days, y = prcp, alpha = 0.001), col = "#73DCFF") +
    ggtitle("Daily precipitation in 2015") +
    xlab("") + ylab("") + 
    theme_light() + scale_alpha(guide = "none")

prcp_2016 <- climate_for_visual %>% filter(year == 2016) %>% 
    select(days, prcp) %>% 
    ggplot(aes(x = days, y = prcp)) + 
    geom_point(col = "blue") + 
    geom_line(aes(x = days, y = prcp, alpha = 0.001), col = "#73DCFF") +
    ggtitle("Daily precipitation in 2016") +
    xlab("") + ylab("") + 
    theme_light() + scale_alpha(guide = "none")

prcp_2017 <- climate_for_visual %>% filter(year == 2017) %>% 
    select(days, prcp) %>% 
    ggplot(aes(x = days, y = prcp)) + 
    geom_point(col = "blue") + 
    geom_line(aes(x = days, y = prcp, alpha = 0.001), col = "#73DCFF") +
    ggtitle("Daily precipitation in 2017") +
    xlab("") + ylab("") + 
    theme_light() + scale_alpha(guide = "none")

# arrange all plots from 2007 to 2012
prcp_07_12 <- ggarrange(prcp_2007, prcp_2008, prcp_2009, 
          prcp_2010, prcp_2011, prcp_2012,
          nrow = 6)
## add same x and y labels
annotate_figure(prcp_07_12, 
                left = text_grob("daily precipitation (mm)", rot = 90),
                bottom = text_grob("consecutive days"))
## publish
#dev.copy(pdf, "prcp0712.pdf")
#dev.off()

# arrange all plots from 2013 to 2018
prcp_13_18 <- ggarrange(prcp_2013, prcp_2014, prcp_2015, 
                        prcp_2016, prcp_2017, prcp_2018,
                        nrow = 6)
## add same x and y labels
annotate_figure(prcp_13_18, 
                left = text_grob("daily precipitation (mm)", rot = 90),
                bottom = text_grob("consecutive days"))
## publish
#dev.copy(pdf, "prcp1318.pdf")
#dev.off()






# Visualize output --------------------------------------------------------
## there are 432 soil loss values, for 36 hillslopes and 12 years.
## plot two hillslopes: Basswood1 hill1 and Orbweaver2 hill3
# ggplots
soil1 <- soil_loss %>% filter(watershed == "Basswood1",
                              hill == "hill1") %>% 
    select(soil_loss, year) %>% 
    ggplot(aes(x = year, y = soil_loss)) + 
    geom_bar(stat = "identity") + 
    ggtitle("Soil loss for Basswood1 hill1") +
    xlab("year") + ylab("soil loss (kg/m^2)") + 
    theme_light()
soil1

soil2 <- soil_loss %>% filter(watershed == "Orbweaver2",
                              hill == "hill3") %>% 
    select(soil_loss, year) %>% 
    ggplot(aes(x = year, y = soil_loss)) + 
    geom_bar(stat = "identity") + 
    ggtitle("Soil loss for Orbweaver2 hill3") +
    xlab("year") + ylab("soil loss (kg/m^2)") + 
    theme_light()

# combine into one
ggarrange(soil1, soil2, nrow = 2)

# publish
#dev.copy(pdf,'soil_loss.pdf')
#dev.copy2pdf(out.type = "pdf")
dev.off()


## soil_loss all in all
soil_loss %>% ggplot(aes(y = soil_loss, x = year)) + 
    geom_boxplot() + ggtitle("Annual soil loss amount for 36 hillslopes") +
    xlab("year") + ylab("soil loss (kg/m^2)") + 
    theme_light()

# publish
#dev.copy(pdf,'soil_loss_all.pdf')
#dev.copy2pdf(out.type = "pdf")
dev.off()






# Visualize hillslope by height and slope ---------------------------------
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

# add hillslope height: standardized height
slope1$height <-  15:1 / slope1$slope
slope1$height[1] <- 550
slope2$height <-  15:1/slope2$slope
slope2$height[1] <- 380

# ggplots of slope
bswd1h1 <- slope1 %>% ggplot(aes(x = n, y = slope)) + 
    geom_step() + ggtitle("Slope profile of Basswood1 hill1") +
    xlab("standardized length") + ylab("slope") + 
    theme_light()

orbw2h3 <- slope2 %>% ggplot(aes(x = n, y = slope)) + 
    geom_step() + ggtitle("Slope profile of Orbweaver2 hill3") +
    xlab("standardized length") + ylab("slope") + 
    theme_light()

# ggplots of heights
bswd1h1_h <- slope1 %>% ggplot(aes(x = n, y = height)) + 
    geom_line() + 
    xlab("standardized length") + ylab("height (m)") + 
    theme_light()

orbw2h3_h <- slope2 %>% ggplot(aes(x = n, y = height)) + 
    geom_line() + 
    xlab("standardized length") + ylab("height (m)") + 
    theme_light()


# combine into one
ggarrange(bswd1h1, bswd1h1_h, nrow = 2, align = "v")
ggarrange(orbw2h3, orbw2h3_h, nrow = 2, align = "v")

# publish
dev.copy(pdf,'bsw_slope.pdf')
dev.off()

dev.copy(pdf,'orb_slope.pdf')
dev.off()



### ADD aLL slopes!!!!
slope_profile_maker <- function(n, d) {
    slopes <- d[n, -c(1, 2)] %>% t() %>% 
        as.data.frame() %>% 
        mutate(nend = c(2:15, NA),
               n = 1:15)
    
    # set correct column names
    colnames(slopes) <- c("slope", "nend", "n") 
    
    # add hillslope height: standardized height
    slopes$height <-  15:1 / slopes$slope
    slopes$height[1] <- slopes$height[2]
    
    # add y end
    slopes$yend <- slopes$slope
    
    plot1 <- slopes %>% 
        ggplot(aes(x = n, y = slope, xend = nend, yend = yend)) +
        geom_vline(aes(xintercept = n), linetype = 2, color = "grey") +
        geom_point() + xlab("") +
        geom_point(aes(x = nend, y = slope), shape = 1) +  
        ggtitle("Slope profile of Orbweaver 3 hill 3") +
        geom_segment() + theme_light()
    
    
    
    plot2 <- slopes %>% ggplot(aes(x = n, y = height)) + 
        geom_line() + 
        geom_vline(aes(xintercept = n), linetype = 2, color = "grey") +
        xlab("standardized length") + ylab("height (m)") + 
        theme_light() 
    
    plot <- ggarrange(plot1, plot2, nrow = 2, align = "v")
    
    return(plot)
}

orb3h3 <- slope_profile_maker(36, num_slope)
orb3h3

ggarrange(orb3h1, orb3h2, orb3h3, nrow = 3)


# publish
dev.copy(pdf, "orb3.pdf")
dev.off()





















