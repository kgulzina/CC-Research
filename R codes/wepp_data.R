######### This code contains WEPP data related information #########
## Author: stat_cat
## Date: 02/04/19





# Sourcing and libraries --------------------------------------------------
source("supplementary_functions.R")
library(devtools)
library(dplyr)
library(stringr)





# STRIPS-1 Configurations -------------------------------------------------
install_github("ISU-STRIPS/STRIPS")
STRIPSMeta::watersheds -> watersheds

### comments: have information on the name and slope of watershed
watersheds





# WEPP watershed visualization --------------------------------------------
# we need to collect data on each watershed/hillslope frame
# frame is given as a total geometric measurements of each hillslope
get_hillslope_frame <- function() {
### write documentation here ###
    # declare names of watersheds
    watersheds <- c("bw1", "bw2", "bw3", "bw4", "bw5", "bw6",
                    "int1", "int2", "int3", "orb1", "orb2", "orb3")
    watersheds_long <- c("Basswood1", "Basswood2", "Basswood3",
                         "Basswood4", "Basswood5", "Basswood6",
                         "Interim1", "Interim2", "Interim3",
                         "Orbweaver1", "Orbweaver2", "Orbweaver3")
    
    # declare hills
    hills <- c("h1", "h2", "h3", "h4")
    hills_long <- c("hill1", "hill2", "hill3", "hill4")
    
    # declare width of each hill
    width <- c(48.0000, 53.0000, 38.0000, 52.0000, 93.0000, 58.0000,
               106.0000, 117.0000, 66.0000, 82.0000, 126.0000, 49.0000)
    
    # declare length of each hill
    length <- c(47.1850, 41.2830, 27.5220, 41.7930, 23.8910, 44.3160,
                14.0000, 51.0880, 29.9120, 23.8820, 59.6930, 35.8350,
                29.0000, 35.1430, 79.0970, 23.8540, 47.6330, 65.5680,
                110.8590, 124.6990, 55.1820, 11.0000, 69.9110, 109.0000,
                81.0890, 29.5590, 94.1690, 5.9800, 125.5980, 37.8720,
                45.5810, 54.0000, 101.1290, 123.0000, 46.0000, 57.0000)
    
    # declare soil type of each hill
    soil_name <- c("LADOGA(SIL)", "LADOGA(SIL)", "LAMONI(L)",
                 "LADOGA(SIL)", "LADOGA(SIL)", "LAMONI(L)",
                 "LADOGA(SIL)", "LADOGA(SIL)", "LAMONI(L)",
                 "LADOGA(SIL)", "LADOGA(SIL)", "LAMONI(L)",
                 "LADOGA(SIL)", "LADOGA(SIL)", "GARA(L)",
                 "LADOGA(SIL)", "LADOGA(SIL)","ARMSTRONG(L)",
                 "OTLEY(SICL)", "CLARINDA(SICL)", "ACKMORE(SIL)",
                 "OTLEY(SICL)", "LADOGA(SIL)", "CLARINDA(SICL)",
                 "ACKMORE(SIL)", "LADOGA(SIL)", "SHELBY(L)",
                 "OTLEY(SICL)", "LADOGA(SIL)", "GARA(L)",
                 "OTLEY(SICL)", "LADOGA(SIL)", "ACKMORE(SIL)",
                 "OTLEY(SICL)", "OTLEY(SICL)", "LADOGA(SIL)")
    soil_id <- c("76C2", "76D2", "822D2", "76C2", "76D2", "822D2",
                 "76C2", "76D2", "822D2", "76C2", "76D2", "822D2",
                 "76C2", "76D2", "993D2", "76C2", "76D2", "993D2",
                 "281C2", "222D2", "5B", "281B", "76C2", "222D2",
                 "5B", "76C2", "93D2", "281C2", "76D2", "179D2",
                 "281C2", "76D2", "5B", "281B", "281C2", "76D2")
    
    # make a count for length and soil type
    count = 1
    
    # dataset frame
    df <- data.frame(HUC_12 = NA, watershed = NA, hill = NA, width = NA,
                     length = NA, soil_name = NA, soil_id = NA,
                     sublength = NA, slope = NA)
    
    for (i in 1:12) {
        # check if interim 2
        if (i == 8) {
            for (j in 1:4) {
                path <- paste("data/slope/071000081505_",
                              watersheds[i], "_", hills[j],
                              ".slp", sep = "")
                print(path)
                d <- read.table(path, skip = 10, header = FALSE)
                
                # length and temp_df are needed for repetitive manipulations
                seq <- seq(1, length(d), 2)
                temp_df <- data.frame(sublength = NA, slope = NA)
                temp_count = 1
                
                # extract slope and sublength of hillslope
                for (k in seq) {
                    temp_df[temp_count,1] <- as.numeric(substr(d[1,k], 1, str_length(d[1,k])-1)) * length[count]
                    temp_df[temp_count,2] <- d[1,k+1] * 100
                    temp_count = temp_count + 1
                }
                
                temp_df <- temp_df %>% mutate(HUC_12 = 071000081505,
                                        watershed = watersheds_long[i],
                                        hill = hills_long[j],
                                        width = width[i],
                                        length = length[count],
                                        soil_name = soil_name[count],
                                        soil_id = soil_id[count])
                # add to the main dataset
                df <- rbind(df, temp_df)
                
                # update count
                count <- count + 1
                
            }
            # check if interim 3
        } else if (i == 9) {
            for (j in 1:2) {
                path <- paste("data/slope/071000081505_",
                              watersheds[i], "_", hills[j],
                              ".slp", sep = "")
                print(path)
                d <- read.table(path, skip = 10, header = FALSE)
                
                # length and temp_df are needed for repetitive manipulations
                seq <- seq(1, length(d), 2)
                temp_df <- data.frame(sublength = NA, slope = NA)
                temp_count = 1
                
                # extract slope and sublength of hillslope
                for (k in seq) {
                    temp_df[temp_count,1] <- as.numeric(substr(d[1,k], 1, str_length(d[1,k])-1)) * length[count]
                    temp_df[temp_count,2] <- d[1,k+1] * 100
                    temp_count = temp_count + 1
                }
                temp_df <- temp_df %>% mutate(HUC_12 = 071000081505,
                                        watershed = watersheds_long[i],
                                        hill = hills_long[j],
                                        width = width[i],
                                        length = length[count],
                                        soil_name = soil_name[count],
                                        soil_id = soil_id[count])
                # add to the main dataset
                df <- rbind(df, temp_df)
                
                # update count
                count <- count + 1
            }
        } else {
            for (j in 1:3) {
                path <- paste("data/slope/071000081505_",
                              watersheds[i], "_", hills[j],
                              ".slp", sep = "")
                print(path)
                d <- read.table(path, skip = 10, header = FALSE)
                
                # length and temp_df are needed for repetitive manipulations
                seq <- seq(1, length(d), 2)
                temp_df <- data.frame(sublength = NA, slope = NA)
                temp_count = 1
                
                # extract slope and sublength of hillslope
                for (k in seq) {
                    temp_df[temp_count,1] <- as.numeric(substr(d[1,k], 1, str_length(d[1,k])-1)) * length[count]
                    temp_df[temp_count,2] <- d[1,k+1] * 100
                    temp_count = temp_count + 1
                }
                
                temp_df <- temp_df %>% mutate(HUC_12 = 071000081505,
                                        watershed = watersheds_long[i],
                                        hill = hills_long[j],
                                        width = width[i],
                                        length = length[count],
                                        soil_name = soil_name[count],
                                        soil_id = soil_id[count])
                # add to the main dataset
                df <- rbind(df, temp_df)
                
                # update count
                count <- count + 1
            }
        }
    }
    return(df[-1,])
}

hillslope_info <- get_hillslope_frame()
hillslope_info %>% glimpse()
hillslope_info %>% head()






# WEPP - climate data - CLIGEN --------------------------------------------
# This climate file was simulated from CLIGEN - Des Moines IA station for 12
# years. I choose Des Moines station, since it is the closest one to the 
# Neal Smith.
climate <- read.table("data/climate/092.63x040.90.cli.txt",
                      skip = 15,
                      header = FALSE,
                      quote = " ")
# although, we have read data correctly there is problem in file format and
# there are only 10 columns, 3 columns are lost. However, we will do not care
# about them, we need only "day", "month", "year" and "precip"
climate %>% head(20)

# dimension is far away from the real dim: 365*12 + 3 = 4383
dim(climate)

# give the column names
colnames(climate) <- c("day",  "month", "year", "prcp", "dur",
                       "tp", "ip", "tmax", "tmin", "rad",
                       "w-vl", "w-dir", "tdew")

### use these lines if you are reading NS real climate file
# write a code which gets rid of unused stuff - 3 columns
# rows to skip:
#skip <- is.na(climate$month)

# I could try to add those three columns
#mess <- climate[skip,1]
#length(mess) 
## 7692 rows which is not the sum of 4383, so I'll not try to fix it!

# valid climate
#climate <- climate[!skip,]

# now, it is correct
#climate %>% head(20)
#dim(climate)

# how many years simulated:
climate$year %>% unique() #12 years -- good






# WEPP - soil input -------------------------------------------------------
## soil encoding -- no header
mapunit <- read.table("data/soil/mapunit.txt",
                      sep = "|",
                      header = FALSE)
mapunit %>% head()

## soil encoding -- no header
muaggatt <- read.table("data/soil/muaggatt.txt",
                       sep = "|")
muaggatt %>% head()

## soil encoding -- no header
muareao <- read.table("data/soil/muareao.txt",
                      sep = "|")
muareao %>% head()

## we will match soil types with id's given here:
muaggatt %>% select(V1, V2, V40) -> soil_id 

## required soil names we see in STRIPS-1 hillslopes
soil_names <- c("76C2", "76D2", "822D2", "993D2", "281C2", "222D2", "5B",
                "281B", "93D2", "179D2")
soil_types <- soil_id[soil_id$V1 %in% soil_names,]
colnames(soil_types) <- c("ID", "name", "MUKEY")
soil_types

# write it to .csv
write.csv(soil_types, "soil_types.csv", row.names = FALSE)

### We have more concrete dataset on soil_types from Gabe:
soil_types_detailed <- read.csv("data/soil/strips_soils_data.txt",
                                header = TRUE,
                                sep = "\t")
soil_types_detailed






# WEPP - soil loss  --------------------------------------------------------
# function to read in
get_env_results <- function() {
### write documentation here ###
    # declare names of watersheds
    watersheds <- c("bw1", "bw2", "bw3", "bw4", "bw5", "bw6",
                    "int1", "int2", "int3", "orb1", "orb2", "orb3")
    watersheds_long <- c("Basswood1", "Basswood2", "Basswood3",
                         "Basswood4", "Basswood5", "Basswood6",
                         "Interim1", "Interim2", "Interim3",
                         "Orbweaver1", "Orbweaver2", "Orbweaver3")
    
    # declare hills
    hills <- c("h1", "h2", "h3", "h4")
    hills_long <- c("hill1", "hill2", "hill3", "hill4")
    
    # dataset frame
    df <- data.frame(HUC_12 = NA, watershed = NA, hill = NA, day = NA, 
                     month = NA, year = NA, Precp = NA, Runoff = NA,
                     IR.det = NA, Av.det = NA, Mx.det = NA, Point = NA, 
                     Av.dep = NA, Max.dep = NA, Point.1 = NA, Sed.Del = NA,
                     ER = NA)
    
    # for faster and correct reading of tables
    col_types <- c("integer", "integer", "integer", "numeric", "numeric", "numeric",
                   "numeric", "numeric", "numeric","numeric", "numeric",
                   "numeric", "numeric", "numeric")
    
    # column names
    col_names <- c("day", "month", "year", "Precp", "Runoff", "IR.det",
                   "Av.det", "Mx.det", "Point", "Av.dep", "Max.dep",
                   "Point.1","Sed.Del", "ER")
    
    for (i in 1:12) {
        # check if interim 1
        if (i == 8) {
            for (j in 1:4) {
                path <- paste("data/env_results/071000081505_",
                              watersheds[i], "_", hills[j],
                              ".env", sep = "")
                print(path)
                data <- read.table(path, skip = 3, header = FALSE,
                                   col.names = col_names,
                                   colClasses = col_types)
                data <- data %>% mutate(HUC_12 = 071000081505,
                                             watershed = watersheds_long[i],
                                             hill = hills_long[j])
                # add to the main dataset
                df <- rbind(df, data)
                
            }
        # check if interim 3
        } else if (i == 9) {
            for (j in 1:2) {
                path <- paste("data/env_results/071000081505_",
                              watersheds[i], "_", hills[j],
                              ".env", sep = "")
                print(path)
                data <- read.table(path, skip = 3, header = FALSE,
                                   col.names = col_names,
                                   colClasses = col_types)
                data <- data %>% mutate(HUC_12 = 071000081505,
                                             watershed = watersheds_long[i],
                                             hill = hills_long[j])
                # add to the main dataset
                df <- rbind(df, data)
            }
        } else {
            for (j in 1:3) {
                path <- paste("data/env_results/071000081505_",
                              watersheds[i], "_", hills[j],
                              ".env", sep = "")
                print(path)
                data <- read.table(path, skip = 3, header = FALSE,
                                   col.names = col_names,
                                   colClasses = col_types)
                data <- data %>% mutate(HUC_12 = 071000081505,
                                             watershed = watersheds_long[i],
                                             hill = hills_long[j])
                # add to the main dataset
                df <- rbind(df, data)
            }
        }
    }
    return(df[-1,])
}

# get the whole dataset
soil_loss_env <- get_env_results()
soil_loss_env %>% glimpse()

# replace years with real ones: 2007-2018
years <- 2007:2018
for (i in 1:12) {
    soil_loss_env$year[soil_loss_env$year == i] = years[i]
}
soil_loss_env %>% glimpse()
soil_loss_env %>% head()

### comments:
# This data contains daily precipitation, sediment loss for years 2007 - 2018
# Days in which there is no sediment loss observed are missing. We can pull
# full daily precipitation data from climate data (above), if needed.





# Combine climate and soil loss datasets ----------------------------------
# get soil loss response as total scalar output
soil_loss_env %>% 
    select(year, watershed, hill, IR.det) %>%
    group_by(year, watershed, hill) %>% 
    summarise(soil_loss = sum(IR.det)) %>% 
    unlist() %>% 
    matrix(ncol = 4, byrow= FALSE) %>% 
    data.frame() -> soil_loss 

# rename columns
colnames(soil_loss) <- c("year", "watershed", "hill", "soil_loss")

# convert soil_loss to numeric
soil_loss %>% 
    mutate(soil_loss = as.numeric(as.character(soil_loss))) -> soil_loss
soil_loss %>% head()
soil_loss %>% glimpse()

### comments: before we get daily precipitation: we have to take into 
# consideration three leap years: 2008, 2012, 2016 which have one additional 
# day. I am removing those days for consistency, if there were no significant 
# precipitation on that day.

# additional days are probably in February, day = 29
leap <- which(climate$month == 2 & climate$day == 29)
climate[leap,] #precipitaion in those days maybe important.. Prof. Niemi?


# get precipitation out of climate
get_prcp <- function(df) {
### write documentation here ###   
    d <- data.frame(t(1:365))
    
    for (i in 2007:2018) {
        temp <- df[df[,1] == i, 2]
        d <- rbind(d, temp)
    }
    
    d <- d[-1,]
    return(d)
}

# add years columns to daily prcp
prcp <- climate[-leap, 3:4] %>% 
    get_prcp() %>% 
    mutate(year = 2007:2018)
prcp %>% head()

### comments:cbind daily prcp and soil loss for each hillslope into one
# we should have three additional columns: field, hill, and soil_loss.
# We will treat 432 examples of soil loss for 12 years and 36 hillslopes as
# the same >> same parameters

# standardize prcp: by taking log(prcp + smallest non-zero prcp)
for (i in 1:12) {
    prcp[i,-366] <- log(prcp[i,-366] + min(prcp[i,prcp[i,-366] != 0]))
}

# we have to merge the datasets by years:
annual_soil_loss <- merge(prcp, soil_loss, by = "year")
annual_soil_loss %>% glimpse()
annual_soil_loss %>% dim()
### comments: I think, I got finally, what I wanted!!!

# save as csv file: annual_soil_loss.csv
write.csv(annual_soil_loss, file = "annual_soil_loss.csv", row.names = FALSE)






# Scalar wepp_data  -------------------------------------------------------
# It's needed to run the simplest model and see if we can get valid
# results. It contains standardized and averaged inputs: temp, prcp, slope
# write a function to convert slope into functional input given: n
gen_slope <- function(n, hillslope_info) {
    # get unique hillslope lengths
    lengths <- unique(hillslope_info$length)
    
    # frame for data
    d <- data.frame(t(1:(n+1))
    for (i in 1:36) {
        # interval length
        h = lengths[i] / n
        
        # which rows to select
        rows <- which(hillslope_info$length == lengths[i])
        
        # how many original sections hillslope has
        M = length(rows)
        
        d[i,0] = 0 #starts from zero
        temp = 0
        for (j in 0:n) {
            temp = temp + h
            for (m in 1:M) {
                if (temp < hillslope_info$sublength[rows[m]]) {
                    ### continue your thing here!!!
                }
            }
        }
    }
}





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










