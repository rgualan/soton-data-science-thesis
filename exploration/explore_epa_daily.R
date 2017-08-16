## Clean environment
rm(list=ls())

## Load libraries
library(openair)

# Read the file ###########################################
d <- read.csv("data/epa/epa_daily/2016/daily_81102_2016_pm10.csv")
d$date <- as.POSIXct(d$Date.Local,format="%Y-%m-%d",tz = "GMT")
d$site <- as.factor(
  paste(d$State.Code,d$County.Code,d$Site.Num,d$POC,substring(d$Local.Site.Name,1,5), sep="-"))
d$var <- d$Arithmetic.Mean

# Filter one state #########################################
d <- d[d$State.Name=="California", ]

sites <- sort(unique(d$site))
index <- 1 
while(index<length(sites)){
  chosen <- sites[index:(index+5)]
  index <- index + 6
  d2 <- d[d$site %in% chosen,]
  summaryPlot(d2, period="months", pollutant = "var", type="site")
  #timePlot(d2, pollutant = "var", type="site")
}


