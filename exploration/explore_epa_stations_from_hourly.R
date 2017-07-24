# Clean environment #################################################################
rm(list=ls())

# Libraries
library(openair)
library(maps) # Provides functions that let us plot the maps
library(mapdata) # Contains the hi-resolution points that mark out the countries.
library(sp)
library(rgeos)
library(ggmap)
library(rgdal)


## Stations monitoring PM10 ####################################################
## Build the list of active stations starting from the data
## Relevant period: 
## 2015-01-01 to 2015-04-30
## Read data
d <- read.csv("data/epa/epa_hourly/2016/hourly_81102_2016_pm10.csv", header=T)
d$date <- as.POSIXct(paste(d$Date.Local,d$Time.Local),format="%Y-%m-%d %H:%M", tz="GMT")
d$site <- as.factor(
  paste(gsub(' ', '.', d$State.Name),d$County.Code,d$Site.Num,d$POC, sep="-"))
d$Station.Code <- as.factor(
  paste(gsub(' ', '.', d$State.Name),d$County.Code,d$Site.Num, sep="-"))
#View(d)
length(unique(d$Station.Code))



## Filter data to match relevant criteria
dateA <- as.POSIXct("2016-01-01", format="%Y-%m-%d", tz="GMT")
dateB <- as.POSIXct("2016-04-30", format="%Y-%m-%d", tz="GMT")
td <- (dateB-dateA)+1
units(td) <- "hours"
td
d2 <- d[d$date>=dateA & d$date<=dateB,]
#View(d2)
#dim(d); dim(d2)


## POCs?
#pocsByStation <- aggregate(POC~site,data=unique(d2$Station.Code,d2$POC),
#                           FUN=length)
#names(pocsByStation)[2] <- "Count"
#View(pocsByStation)
# It gets stalled for some reason!


## Number of records by station
recordsByStation <- aggregate(Sample.Measurement~Station.Code,data=d2,FUN=length)
names(recordsByStation)[2] <- "Count"
#View(recordsByStation)

#Note:
#The following stations seem to have two POCS:
#Wisconsin-111-7 California-19-11 Iowa-163-17 Missouri-97-3 New.Mexico-13-20 


## Plot some sites (Checking NAs)
summaryPlot(d2[d2$Station.Code %in% c("Illinois-31-1016","Nevada-3-1019",
                              "California-19-500","New.Mexico-13-20"),
               c("date","site","Sample.Measurement")],
            period = "months", pollutant="Sample.Measurement", type="site")
readline("Continue?")


## Active stations
minNumRecords <- 0.95*as.integer(td)
stations <- recordsByStation[recordsByStation$Count>minNumRecords,]
#View(stations)
## Complete fields of the stations
stations_md <- read.csv("data/epa/sites/aqs_sites.csv", header=T)
stations_md <- unique(stations_md[,c("State.Code","County.Code",
                                     "Site.Number","Latitude","Longitude",
                                     "Datum","Elevation", "Land.Use",
                                     "Location.Setting", "State.Name",
                                     "County.Name")])
stations_md$Station.Code <- as.factor(
  paste(gsub(' ', '.', stations_md$State.Name),stations_md$County.Code,
        stations_md$Site.Num, sep="-"))

stations <- merge(stations,stations_md,all.x=T)
#View(stations)
nrow(stations)



# Plot all the active stations in USA
spStations <- stations
coordinates(spStations) <- ~Longitude+Latitude
proj4string(spStations) <- "+proj=longlat +datum=WGS84"
map('state')	
points(spStations, cex=0.5)
title("Active hourly PM10 stations")
readline("Continue?")

#mapUSA <- get_map(location = 'USA', zoom = 4, maptype = "roadmap")
#save(mapUSA, file="data/maps/mapUSA.Rdata")
if(!exists("mapUSA")) load("data/maps/mapUSA.Rdata")

ggmap(mapUSA) +
  ggtitle("Active hourly PM10 stations") +
  geom_point(aes(x=Longitude, y=Latitude, 
                 fill = Land.Use), 
             data = stations, alpha = .75, shape=24, size=1)
nrow(stations)


# Focus on California
stations2 <- stations[stations$State.Name=="California",]
map('state', 'California')	
points(Latitude~Longitude,stations2,cex=2,pch=".")
title("Active hourly PM stations in California")
readline("Continue?")

# Stations by Land Use
#mapCali <- get_map(location = 'California', zoom = 6, maptype = "roadmap")
#save(mapCali, file="data/maps/mapCali.Rdata")
if(!exists("mapCali")) load("data/maps/mapCali.Rdata")
ggmap(mapCali) +
  ggtitle("Active hourly PM10 stations in California") +
  geom_point(aes(x=Longitude, y=Latitude, 
                 fill = Land.Use), 
             data = stations2, alpha = .75, shape=24, size=1)
readline("Continue?")

# Stations by Location Setting
ggmap(mapCali) +
  ggtitle("Active hourly PM10 stations in California") +
  geom_point(aes(x=Longitude, y=Latitude, 
                 fill = Location.Setting), 
             data = stations2, alpha = .75, shape=24, size=1)
readline("Continue?")





