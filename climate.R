### observed climate data
## Keokuk/Iowa/Daily

keukok <- read.csv(file = "Keokuk Iowa NOAA Daily.txt")
dim(keukok)
summary(keukok)
# some measurements in precipitation is recorded as logical??

unique(keukok$year) # 13 years
tapply(as.numeric(keukok$Prcp), INDEX = keukok$year, mean)
# highest prcp is for 1996: some anomaly?








































