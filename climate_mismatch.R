######### This code contains WEPP data related information #########
## Author: stat_cat
## Date: 05/21/19





# Sourcing and libraries --------------------------------------------------
library(devtools)
library(dplyr)
library(stringr)





# Handling mismatch example -----------------------------------------------
climate_daily <- read.table("data/climate_mismatch/weppdemo_daily.cli",
                            skip = 15,
                            header = FALSE,
                            quote = " ")
climate_daily %>% head(20)

# dimension is far away from the real dim: 365*12 + 3 = 4383
dim(climate_daily)

# give the column names
colnames(climate_daily) <- c("day",  "month", "year", "prcp", "dur",
                             "tp", "ip", "tmax", "tmin", "rad", "w-vl", "w-dir", "tdew")

# select the first year only
climate_daily <- climate_daily %>% 
    select(day, month, year, prcp) %>% 
    filter(year == 2007)


# read env_result
soil_daily <- read.table("data/climate_mismatch/daily.txt",
                         skip = 3)
colnames(soil_daily) <- c("day", "month", "year", "Precp", "Runoff", "IR.det",
                          "Av.det", "Mx.det", "Point", "Av.dep", "Max.dep",
                          "Point.1","Sed.Del", "ER")

# select only needed part and rename year
soil_daily <- soil_daily %>% 
    select(day, month, year, Precp, Runoff, IR.det) %>% 
    mutate(year = 2007)
soil_daily

# merge into one
true_match <- merge(climate_daily, soil_daily, by = c("day", "month", "year"))
true_match

# write into csv
write.csv(true_match, file = "truematch.csv")





# Comments: ---------------------------------------------------------------
# I looked to the wepp usernum more closely and saw that there are some
# differences between NS.cli.txt and the default climate file generated 
# by CLIGEN.

# Second line in climate txt input file contains three integers:
# a) simulation mode - integer (itemp)
#   1 - continuous
#   2 - single storm
# b) breakpoint data flag - integer (ibrkpt)
#   0 - no breakpoint data used
#   1 - breakpoint data used
# c) wind information/ET equation flag - integer (iwind)
#   0 - wind information exists - use Penman ET equation
#   1 - no wind information exists - use Priestley-Taylor ET equation

## cligen 2nd line      = 1 0 0
## NC.cli.txt 2nd line  = 1 1 0
## So the differenceas are because of breakpoint data used in NS.cli.txt

# More detailed information can be found on page 10 of wepp_usernum.pdf