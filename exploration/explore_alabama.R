# Clean environment
rm(list=ls())

# Libraries
library(openair)

# Read the data
df <- read.csv('data/epa_hourly/hourly_81102_2017_pm10_alabama.csv',header=T)
df$DateTime.Local <- paste(df$Date.Local,df$Time.Local)
df$DateTime.Local <- as.POSIXct(df$DateTime.Local,format="%Y-%m-%d %H:%M",tz="GMT")
summary(df)

# Simple plot
df1 <- df[df$Site.Num==1,c('DateTime.Local','Sample.Measurement','Site.Num')]
summary(df1)
plot(df1$DateTime.Local, df1$Sample.Measurement, type='l')


# Multi-site plot
#df2 <- df[df$Site.Num==1,c('DateTime.Local','Sample.Measurement','Site.Num')]
df2 <- df[df$Site.Num>=1 & df$Site.Num<=10,
          c('DateTime.Local','Sample.Measurement','Site.Num')]
df2$date <- df2$DateTime.Local
df2$pm10 <- df2$Sample.Measurement
df2$site <- as.factor(paste('Site-',df2$Site.Num))
str(df2)
a <- df2[,c('site','date','pm10')]
summaryPlot(a, type='site', pollutant = "pm10", period = "months")


