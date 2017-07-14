library(openair)

# Get the file #############################################
# download.file("https://aqsdr1.epa.gov/aqsweb/aqstmp/airdata/daily_44201_2017.zip",
#               destfile="data/epa_hourly/daily.zip")
# unzip("data/epa_hourly/daily.zip", exdir="data/epa_hourly")

# Check the file ###########################################
d <- read.csv("data/epa_hourly/daily_44201_2017.csv")

# Filter one state #########################################
d <- d[d$State.Name=="Utah", ]
d$date <- as.POSIXct(d$Date.Local,format="%Y-%m-%d",tz = "GMT")
d$site <- as.factor(
  paste(d$State.Code,d$County.Code,d$Site.Num,d$POC,substring(d$Local.Site.Name,1,5), sep="-"))
d$ozone <- d$Arithmetic.Mean

summaryPlot(d, period="months", pollutant = "ozone", type="site")
readline("Continue?")
timePlot(d, pollutant = "ozone", type="site")
