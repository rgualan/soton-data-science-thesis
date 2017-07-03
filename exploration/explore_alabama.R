library(openair)

# Read the data
df <- read.csv('data/epa_hourly/hourly_81102_2017_pm10_alabama.csv',header=T)
df$DateTime.Local <- paste(df$Date.Local,df$Time.Local)
df$DateTime.Local <- strptime(df$DateTime.Local,format="%Y-%m-%d %H:%M")
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
summaryPlot(a, type='site')


# mydata
mydata <- read.csv("data/openair/OpenAir_example_data_long.csv", header = TRUE)
mydata$date <- as.POSIXct(strptime(mydata$date, format = "%d/%m/%Y %H:%M", tz = "GMT"))
head(mydata)
plot(mydata$date, mydata$pm10, type = 'line')
hist(mydata$no2)
summaryPlot(mydata)
summaryPlot(mydata, clip = FALSE) # do not clip density plot data
summaryPlot(mydata, percentile = 0.95) # exclude highest 5 % of data etc.
summaryPlot(mydata, na.len = 10)
summaryPlot(mydata, col.data = "green") # show data in green
summaryPlot(mydata, col.dens = "black") # show density plot line in black

pairs(mydata[sample(1:nrow(mydata), 500), c(1, 2, 3, 4, 5)],
      lower.panel = panel.smooth,
      upper.panel = NULL,
      col = "skyblue3")

timePlot(mydata, pollutant = c("nox", "no2", "o3", "pm10"))
timePlot(selectByDate(mydata, year = 2003, month = "aug"),
         pollutant = c("nox", "o3", "pm25", "pm10", "ws"))

calendarPlot(mydata, pollutant = "o3", year =2003)

TheilSen(mydata, pollutant = "o3", ylab = "ozone (ppb)", deseason = TRUE)

#par(mfrow=c(2,2)) # It does not work with smoothTrend
smoothTrend(mydata, pollutant = "o3", ylab = "concentration (ppb)",
            main = "monthly mean o3")
smoothTrend(mydata, pollutant = "o3", deseason = TRUE, ylab = "concentration (ppb)",
            main = "monthly mean deseasonalised o3")
smoothTrend(mydata, pollutant = "no2", simulate = TRUE, ylab = "concentration (ppb)",
            main = "monthly mean no2 (bootstrap uncertainties)")
smoothTrend(mydata, pollutant = "no2", deseason = TRUE, simulate =TRUE,
            ylab = "concentration (ppb)",
            main = "monthly mean deseasonalised no2 (bootstrap uncertainties)")

timeVariation(subset(mydata, ws > 3 & wd > 100 & wd < 270),
              pollutant = "pm10", ylab = "pm10 (ug/m3)")

data2003 <- selectByDate(mydata, year = 2003)
scatterPlot(data2003, x = "nox", y = "no2")
scatterPlot(selectByDate(mydata, year = 2003), x = "nox", y = "no2",
            method = "density", col = "jet")
scatterPlot(data2003, x = "nox", y = "no2", z = "o3", type = c("season", "weekend"),
            limits = c(0, 30))



### Google Maps
load(url("http://www.erg.kcl.ac.uk/downloads/Policy_Reports/AQdata/o3Measurements.RData"))
head(o3Measurements)
summary(o3Measurements)
load(url("http://www.erg.kcl.ac.uk/downloads/Policy_Reports/AQdata/siteDetails.RData"))
head(siteDetails)

## cut data into seasons
## load plyr package
library(plyr)
o3Measurements <- cutData(o3Measurements, "season")
## calculate means/maxes and merge...
annual <- ddply(o3Measurements, .(site), numcolwise(mean), na.rm = TRUE)
## by site AND season
means <- ddply(o3Measurements, .(site, season), numcolwise(mean), na.rm = TRUE)
peaks <- ddply(o3Measurements, .(site, season), numcolwise(max), na.rm = TRUE)
annual <- merge(annual, siteDetails, by = "site")
means <- merge(means, siteDetails, by = "site")
peaks <- merge(peaks, siteDetails, by = "site")
## now make first plot
googleMapsPlot(annual, lat = "latitude", long = "longitude", pollutant = "o3",
               maptype = "roadmap", col = "jet")


# AURN dataset
mary <- importAURN(site = "my1", year = 2010:2016)
abd <- importAURN(site="ABD7", year = 2010:2016)
summaryPlot(mary)
summaryPlot(mary, pollutant = c('pm10'))
summaryPlot(abd)

# KCL dataset
a30 <- importAURN(site="A30", year = 2010:2016) # Not working???

# AIR
kc1.airbase <- importAirbase(site = "gb0620a")
summary(kc1.airbase)
head(kc1.airbase)
summaryPlot(kc1.airbase)
