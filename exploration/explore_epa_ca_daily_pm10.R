## Clean environment #################################################################
rm(list=ls())
par(ask=T)

## Libraries
library(lattice)
library(RColorBrewer)
library(ggplot2)

## Read data #########################################################################
## Measurementiable: PM10
## Build the list of active stations starting from the data
## Relevant period: 
## Read data
d <- read.csv("data/epa/epa_daily/2016/daily_81102_2016_pm10.csv", header=T)
d$Date <- as.POSIXct(d$Date.Local,format="%Y-%m-%d", tz="GMT")
## Filter California
d <- d[d$State.Name=="California", ]
## Station code
## Option a
# d$Site <- as.factor(
#   paste(gsub(' ', '.', d$State.Name),d$County.Code,d$Site.Num,d$POC, sep="-"))
## Option b
d$Site <- as.factor(sprintf("%03d-%04d-%01d",d$County.Code,d$Site.Num,d$POC))
# Code a)
# d$Station.Code <- as.factor(
#   paste(gsub(' ', '.', d$State.Name),d$County.Code,d$Site.Num, sep="-"))
# Code b) SS-CCC-NNNN-PPPPP-Q
d$Station.Code <- as.factor(sprintf("%03d-%04d",d$County.Code,d$Site.Num))
## Generic variable name
names(d)[names(d)=="Arithmetic.Mean"] <- "Measurement"
#View(d)
## Important fields
d <- d[,c("Station.Code","Site","Latitude","Longitude","Date","Measurement")] # "Date","State.Code","State.Name", "County.Code","Site.Num","POC",
head(d)

## Apply filters 
dateA <- as.POSIXct("2016-01-01", format="%Y-%m-%d", tz="GMT")
dateB <- as.POSIXct("2016-12-31", format="%Y-%m-%d", tz="GMT")
td <- (dateB-dateA)+1
# d <- d[d$Date>=dateA & d$Date<=dateB, ]

## Read sites
sites_md <- read.csv("data/epa/sites/aqs_sites.csv", header=T)
sites_md <- sites_md[sites_md$State.Name=="California",]
sites_md$Station.Code <- as.factor(sprintf("%03d-%04d",sites_md$County.Code,sites_md$Site.Num))
#View(sites_md)
## Combine sites (data) with sites (metadata)
sites <- sort(unique(d$Station.Code))
sitesDs <- merge(data.frame(Station.Code=sites), 
                 sites_md[,c("Station.Code","Latitude","Longitude","Datum","Elevation",
                             "Location.Setting")])
#View(sitesDs)

## Filter only rural stations
d <- merge(d, sitesDs[,c("Station.Code","Location.Setting")])
head(d)
unique(d$Location.Setting)
# d <- d[d$Location.Setting!="URBAN AND CENTER CITY",]
# unique(d$Station.Code)

## Pre-process ##############################################################################
## Check POC 
# tbl <- table(d[,c("Station.Code","POC")])
# tbl <- cbind(tbl,rowSums(tbl))
# tbl[tbl[,8]>0,]
# ## Check Sample.duruation
# tbl <- table(d[,c("Station.Code","Sample.Duration")])
# tbl[tbl[,1]>0 | tbl[,2]>0,]
# ## Check DATUM
# tbl <- table(d[,c("Station.Code","Datum")])
# tbl[tbl[,1]>0 | tbl[,2]>0,] ## All WGS84
# ## Check Event.Type
# tbl <- table(d[,c("Station.Code","Event.Type")])
# tbl[rowSums(tbl)>0,]
# d <- d[d$Event.Type != "Included",]
## Remove erratic station
# d[d$Measurement == max(d$Measurement),]
# d <- d[!d$Site %in% c("California-51-11-3","California-27-25-2","California-27-29-1",
#                       "California-65-2005-3", "California-25-7-3"),]
maxBySite <- aggregate(Measurement~Site,d,max)
maxBySite <- maxBySite[order(maxBySite$Measurement, decreasing = T),]
d <- d[!d$Site %in% head(maxBySite$Site, n=16),]
ggplot(d) + 
  geom_line(aes(x=Date,y=Measurement,col=Station.Code), alpha=0.5, size=0.5) +
  theme(legend.position="none")
## Negative values?
d[d$Measurement<0,]
d[d$Measurement<0,]$Measurement<-0
## Summary
hist(d$Measurement)
range(d$Measurement)


## Melt POC
dim(d)
d <- aggregate(Measurement~Station.Code+Date,d,mean)
dim(d)

## Check NAs ##########################################################
## Plot a data concentration as matrix (sites x date)
NS<-50
for(i in 1:ceiling(length(sites)/50)){
  #i = 1
  print(i)
  print(levelplot(Measurement~Date*Station.Code, 
                  d[d$Station.Code %in% sites[(NS*(i-1)+1):(NS*i)],],
                  cuts=10,col.regions=rev(brewer.pal(11,"Spectral")),
                  scales=list(y=list(cex=.7))))
  
}

## Disable stations with no enough data ################################
recordsByStation <- aggregate(Measurement~Station.Code, d, length)
names(recordsByStation)[ncol(recordsByStation)] <- "Count"
#View(recordsByStation);
minNumRecords <- 0.70*as.integer((dateB-dateA)+1)
## Disable stations
## Before:
nrow(recordsByStation)
## After 
sites.on <- recordsByStation[recordsByStation$Count>minNumRecords,]
nrow(sites.on)
# table(merge(sites.on,sites_md)$Location.Setting)
## Apply filter
d <- d[d$Station.Code %in% sites.on$Station.Code,]
sitesDs <- sitesDs[sitesDs$Station.Code %in% sites.on$Station.Code,]
sitesDs[sitesDs$Location.Setting=="",]$Location.Setting <- "RURAL"

## Plot active stations ##################################################################
## Map of California
mapUSA <- readRDS("data/maps/usa/USA_adm1.rds")
mapCA <- mapUSA[mapUSA$NAME_1=="California",]
proj4string(mapCA)
plot(mapCA)

ggplot(mapCA) +
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", colour = "black") +
  geom_point(data = sitesDs, aes(x = Longitude, y = Latitude, fill=Location.Setting),
             alpha = .75, shape=21, size=2) +
  labs(x = "Longitude", y = "Latitude", title = "Map of California") +
  coord_quickmap()


## Check NAs and Time series ##################################################################
## Concentration map
print(levelplot(Measurement~Date*Station.Code,d,
                cuts=10,col.regions=rev(brewer.pal(11,"Spectral")),
                scales=list(y=list(cex=.7))))


## Time series plot
ggplot(d) + 
  geom_line(aes(x=Date,y=Measurement,col=Station.Code), alpha=0.5, size=0.5) +
  theme(legend.position="none")

## Notes:
## Maybe remove the station with the event! 
range(d$Measurement)
hist(d$Measurement)
max(d$Measurement)

## Save dataset
saveRDS(d, "data/epa/epa_daily/2016/california.RDS")
saveRDS(sitesDs, "data/epa/epa_daily/2016/sites.RDS")


## End
par(ask=F)