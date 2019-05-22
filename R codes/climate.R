### observed climate data
## Keokuk/Iowa/Daily


# Libraries ---------------------------------------------------------------
library(dplyr)
library(lubridate)
library(ggplot2)


keukok <- read.csv(file = "Keokuk Iowa NOAA Daily.txt")
dim(keukok)
summary(keukok)
# some measurements in precipitation is recorded as logical??

unique(keukok$year) # 13 years
tapply(as.numeric(keukok$Prcp), INDEX = keukok$year, mean)
# highest prcp is for 1996: some anomaly?

# mutate new variables: dates and total precipitation
# precipitaions were recorded as factors, as well as there are some
# values as T: maybe true? they couldn't measure the precipitation,
# but it snowed or rained on that day for sure?
# For now, i will just keep those values as 0 >> for future consideration

# did not allow do it by dplyr
keukok$Prcp[which(keukok$Prcp == "T")] = as.factor(0)
keukok$Prcp <- as.numeric(as.character(keukok$Prcp))
keukok$Prcp[which(keukok$Prcp == 999.99)] = 0

# extract needed parts
keukok %>% select(year, month, day, Prcp) -> new_keukok


# some insights
plot(new_keukok$Prcp)
abline()

# try to see the periodicity -- have to rely on this for now
new_keukok %>% filter(year > 1999) %>% ggplot(aes(month, Prcp)) + geom_point() + facet_grid(.~year)





# Climate Data for Iowa ---------------------------------------------------
iowa <- read.csv("climate-trends---state_-ia,-season_-annual.csv")
iowa %>% glimpse()
## yearly data, we need monthly -- even daily









































