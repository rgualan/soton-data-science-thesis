# Clean environment #################################################################
rm(list=ls())

# Libraries #########################################################################
library(openair)


# Read the data #####################################################################
d <- read.csv('data/epa/epa_hourly/2016/hourly_42101_2016_co.csv',header=T)
d$DateTime.Local <- paste(d$Date.Local,d$Time.Local)
d$date <- as.POSIXct(d$DateTime.Local,format="%Y-%m-%d %H:%M", tz="GMT")
d$site <- as.factor(
  paste(d$State.Code,d$County.Code,d$Site.Num,d$POC,substring(d$Local.Site.Name,1,5), sep="-"))
d$var <- d$Sample.Measurement
simpleFilter <- head(unique(d$site))
summaryPlot(d[d$site %in% simpleFilter,c("date","site","var")], type="site")


