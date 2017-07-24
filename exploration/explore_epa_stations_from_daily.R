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


## Stations #########################################################################
## Variable: PM10
## Build the list of active stations starting from the data
## Relevant period: 
## 2015-01-01 to 2015-04-30
## Read data
varName="Pm10"
varName="Ozone"
#d <- read.csv("data/epa/epa_daily/2016/daily_81102_2016_pm10.csv", header=T)
d <- read.csv("data/epa/epa_daily/2016/daily_44201_2016_ozone.csv", header=T)
d$date <- as.POSIXct(d$Date.Local,format="%Y-%m-%d", tz="GMT")
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
recordsByStation <- aggregate(Arithmetic.Mean~Station.Code,data=d2,FUN=length)
names(recordsByStation)[2] <- "Count"
#View(recordsByStation)


## Plot time sereies of random stations (Overall evaluation of NAs)
randomStations <- recordsByStation[sample(nrow(recordsByStation),5),c("Station.Code")]
summaryPlot(d2[d2$Station.Code %in% randomStations,
               c("date","site","Arithmetic.Mean")],
            period = "months", pollutant="Arithmetic.Mean", type="site",
            main="Time series of random stations")
readline("Continue?")


## Active stations
minNumRecords <- 0.90*as.integer(td)
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
map('state')	
points(Latitude~Longitude,stations,cex=0.5,pch=17,col="blue")
title(sprintf("Active daily %s stations",varName))
readline("Continue?")

#mapUSA <- get_map(location = 'USA', zoom = 4, maptype = "roadmap")
#save(mapUSA, file="data/maps/mapUSA.Rdata")
if(!exists("mapUSA")) load("data/maps/mapUSA.Rdata")

ggmap(mapUSA) +
  ggtitle(sprintf("Active daily %s stations",varName)) +
  geom_point(aes(x=Longitude, y=Latitude, fill=Location.Setting), 
             data = stations, alpha = .75, shape=24, size=1)
nrow(stations)
readline("Continue?")

## California
stations2 <- stations[stations$State.Name=="California",]
# map('state', 'California')	
# points(Latitude~Longitude,stations2,cex=2,pch=".")
# title(sprintf("Active hourly %s stations in California",varName))
# readline("Continue?")
# nrow(stations2)

## Map
if(!exists("mapCali")) load("data/maps/mapCali.Rdata")
ggmap(mapCali) +
  ggtitle(sprintf("Active hourly %s stations in California",varName)) +
  geom_point(aes(x=Longitude, y=Latitude, fill=Location.Setting), 
             data = stations2, alpha = .75, shape=24, size=1)
readline("Continue?")

## New York
stations2 <- stations[stations$State.Name=="New York",]
#mapNY <- get_map(location = 'New York', zoom = 6, maptype = "roadmap")
#save(mapNY, file="data/maps/mapNY.Rdata")
if(!exists("mapNY")) load("data/maps/mapNY.Rdata")
ggmap(mapNY) +
  ggtitle(sprintf("Active hourly %s stations in NY",varName)) +
  geom_point(aes(x=Longitude, y=Latitude, fill=Location.Setting), 
             data = stations2, alpha = .75, shape=24, size=2)
nrow(stations2)
readline("Continue?")




## Check Daily PM2.5 ####################################################################
## PM2.5 has two datasets: FRM/FEM and non FRM/FEM
## Which might mean <accurate> and <inaccurate>
## Read data
da <- read.csv("data/epa/epa_daily/2016/daily_88101_2016_pm25_FRMFEM.csv", header=T)
db <- read.csv("data/epa/epa_daily/2016/daily_88502_2016_non_FRMFEM.csv", header=T)
da$FRMFEM <- "FRMFEM"
db$FRMFEM <- "non-FRMFEM"
#View(da);View(db) 
d <- rbind(da,db)

d$date <- as.POSIXct(d$Date.Local,format="%Y-%m-%d", tz="GMT")
d$site <- as.factor(
  paste(gsub(' ', '.', d$State.Name),d$County.Code,d$Site.Num,d$POC, sep="-"))
d$Station.Code <- as.factor(
  paste(gsub(' ', '.', d$State.Name),d$County.Code,d$Site.Num, sep="-"))
#View(d)
length(unique(d$Station.Code))


# Check what months have the most data
dCount <- d
dCount$MonthYear <- as.factor(format(d$date, format="%Y-%m"))
countByMonthYear <- aggregate(Arithmetic.Mean~MonthYear,data=dCount,FUN=length)
names(countByMonthYear)[2] <- "Count"
#View(countByMonthYear)
par(las=2) # make label text perpendicular to axis
barplot(countByMonthYear$Count, names.arg=countByMonthYear$MonthYear, 
        horiz=T, cex.names=0.7, main="Data by month") 
readline("Continue?")


## Plot all the stations in USA
stations <- unique(d[,c("Station.Code","Longitude","Latitude","FRMFEM")])
map('state')	
points(Latitude~Longitude,stations[stations$FRMFEM=="FRMFEM",],cex=0.5,pch=16,col="blue")
points(Latitude~Longitude,stations[stations$FRMFEM=="non-FRMFEM",],cex=0.5,pch=16,col="green")
title("Daily PM2.5 stations FEM versus non FEM")
readline("Continue?")


## Plot ACTIVE stations in USA 
## Filter data to match relevant criteria
d2 <- d[d$date>=dateA & d$date<=dateB,]
#View(d)
## How many records per stations
recordsByStation <- aggregate(Arithmetic.Mean~Station.Code+Longitude+Latitude+FRMFEM, 
                              data=d2,FUN=length)
names(recordsByStation)[5] <- "Count"
#View(recordsByStation);
minNumRecords <- 0.90*as.integer((dateB-dateA)+1)
## Active stations
stations <- recordsByStation[recordsByStation$Count>minNumRecords,]
nrow(stations)
## Plot 
if(!exists("mapUSA")) load("data/maps/mapUSA.Rdata")
ggmap(mapUSA) +
  ggtitle("Active daily PM2.5 stations") +
  geom_point(aes(x=Longitude, y=Latitude, fill=FRMFEM), 
             data=stations, alpha=.75, shape=24, size=1)




## Check Daily Temperature ####################################################################
## Read data
d <- read.csv("data/epa/epa_daily/2016/daily_TEMP_2016.csv", header=T)
d$date <- as.POSIXct(d$Date.Local,format="%Y-%m-%d", tz="GMT")
d$site <- as.factor(
  paste(gsub(' ', '.', d$State.Name),d$County.Code,d$Site.Num,d$POC, sep="-"))
d$Station.Code <- as.factor(
  paste(gsub(' ', '.', d$State.Name),d$County.Code,d$Site.Num, sep="-"))
#View(d)
length(unique(d$Station.Code))



## Filter data to match relevant criteria
d2 <- d[d$date>=dateA & d$date<=dateB,]
#View(d2)
#dim(d); dim(d2)


## Number of records by station
recordsByStation <- aggregate(Arithmetic.Mean~Station.Code,data=d2,FUN=length)
names(recordsByStation)[2] <- "Count"
#View(recordsByStation)


## Active stations
minNumRecords <- 0.95*as.integer(dateB-dateA+1)
stations <- recordsByStation[recordsByStation$Count>minNumRecords,]
#View(stations)
nrow(stations)

## Complete metadata fields
stations <- merge(stations,stations_md,all.x=T)

## Plot the active stations
if(!exists("mapUSA")) load("data/maps/mapUSA.Rdata")

ggmap(mapUSA) +
  ggtitle("Active daily Temperature stations") +
  geom_point(aes(x=Longitude, y=Latitude, fill = Location.Setting), #Location.Setting/Land.Use
             data = stations, alpha = .75, shape=24, size=1)
