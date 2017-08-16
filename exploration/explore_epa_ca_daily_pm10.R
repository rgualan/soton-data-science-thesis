## Clean environment #################################################################
rm(list=ls())

## Libraries
library(lattice)
library(RColorBrewer)
library(ggplot2)
source("util/my_helper.R")

## Global variables
paper <- setupPaper()
numDays <- as.integer((convertStringToPOSIXct("2016-12-31")-convertStringToPOSIXct("2016-01-01"))+1)

## Basic pre-processing #########################################################################
if(F){
  ## Build the list of active stations starting from the data
  # d <- read.csv("data/epa/epa_daily/2016/daily_81102_2016_pm10.csv", header=T)
  # d$Date <- as.POSIXct(d$Date.Local,format="%Y-%m-%d", tz="GMT")
  # ## Filter California
  # d <- d[d$State.Name=="California", ]
  # d$Site <- as.factor(sprintf("%03d-%04d-%01d",d$County.Code,d$Site.Num,d$POC))
  # d$Station.Code <- as.factor(sprintf("%03d-%04d",d$County.Code,d$Site.Num))
  # ## Generic variable name
  # names(d)[names(d)=="Arithmetic.Mean"] <- "Measurement"
  # ## Keep relevant fields
  # d <- d[,c("Station.Code","Site","POC", "Latitude","Longitude","Date","Measurement")] # "Date","State.Code","State.Name", "County.Code","Site.Num","POC",
  # head(d)
  # saveRDS(d, "data/epa/epa_daily/2016/daily_81102_2016_pm10_ca.RDS")
  d <- readRDS("data/epa/epa_daily/2016/daily_81102_2016_pm10_ca.RDS")
  sites <- getSites(d)
  
  ## Filter only rural stations ## WARNING!!!
  ## Are this filters necessary???
  # d <- merge(d, sitesDs[,c("Station.Code","Location.Setting")])
  # head(d)
  # unique(d$Location.Setting)
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
  # d[d$PM10 == max(d$PM10),]
  # d <- d[!d$Site %in% c("California-51-11-3","California-27-25-2","California-27-29-1",
  #                       "California-65-2005-3", "California-25-7-3"),]
  maxBySite <- aggregate(Measurement~Site,d,max)
  maxBySite <- maxBySite[order(maxBySite$Measurement, decreasing = T),]
  d <- d[!d$Site %in% head(maxBySite$Site, n=16),]
  ## Negative values?
  d[d$Measurement<0,]
  d[d$Measurement<0,]$Measurement<-0
  ## Summary

  ## Melt POC
  dim(d)
  d <- aggregate(Measurement~Station.Code+Date,d,mean)
  dim(d)

  ## Check NAs ##########################################################
  ## Plot a data concentration as matrix (sites x date)
  NS<-50
  if(F){
    for(i in 1:ceiling(length(sites)/50)){
      #i = 1
      print(i)
      print(levelplot(Measurement~Date*Station.Code, 
                      d[d$Station.Code %in% sites[(NS*(i-1)+1):(NS*i)],],
                      cuts=10,col.regions=rev(brewer.pal(11,"Spectral")),
                      scales=list(y=list(cex=.7))))
      
    }
  }
  
  ## Disable stations with no enough data ################################
  recordsByStation <- aggregate(Measurement~Station.Code, d, length)
  names(recordsByStation)[ncol(recordsByStation)] <- "Count"
  #View(recordsByStation);
  minNumRecords <- 0.80*numDays
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

  ## Notes:
  ## Maybe remove the station with the event! 
  range(d$Measurement)
  hist(d$Measurement)
  max(d$Measurement)

  ## Save dataset
  names(d)[3] <- "PM10" 
  saveRDS(d, "data/epa/epa_daily/2016/california_pm10.RDS")
}

## Read data ###########################################################################
d <- readRDS("data/epa/epa_daily/2016/california_pm10.RDS")

## Initial exploration ##############################################################################
printPlot(paper,"img/eda/hist_pm10.jpeg",5,5,FUN=function(){
  hist(d$PM10, main="", xlab="PM10")  
})
printPlot(paper,"img/eda/hist_pm10_2.jpeg",5,5,FUN=function(){
  hist(log(d$PM10), main="", xlab="log(PM10)")  
})
range(d$PM10)

if(F){
  ggplot(d) + 
    geom_line(aes(x=Date,y=PM10,col=Station.Code), alpha=0.5, size=0.5) +
    theme(legend.position="none")
}







## Plot active stations ##################################################################
printPlot(paper,"img/eda/stations_ca_pm10.jpeg",6,6,FUN= function(){
  p <- ggplot(getCAmap()) +
    geom_polygon(aes(x = long, y = lat, group = group), fill = "white", colour = "black") +
    geom_point(data = sitesDs, aes(x = Longitude, y = Latitude, fill=Location.Setting),
               alpha = .75, shape=21, size=2) +
    labs(x = "Longitude", y = "Latitude") +
    coord_quickmap()
})

## Check NAs and Time series ##################################################################
## Concentration map
# print(levelplot(PM10~Date*Station.Code,d,
#                 cuts=10,col.regions=rev(brewer.pal(11,"Spectral")),
#                 scales=list(y=list(cex=.7))))
printPlot(paper, "img/eda/heatmap_ca_pm10.jpeg", 6, 7, FUN= function(){
  heatmapPlusTs(d, "PM10", "PM10 (ppb)")
})


## Statistics about NAs
recordsByStation <- aggregate(PM10~Station.Code, d, length)
names(recordsByStation)[2] <- "Count" 
recordsByStation$Missing <- numDays-recordsByStation$Count
recordsByStation$MissingP <- (recordsByStation$Missing/numDays)*100
mean(recordsByStation$MissingP)
range(recordsByStation$MissingP)

## Statistics about Distance
sites <- getSites(d)
sites.sp <- convertDataToSp(sites)
sites.d <- dist(coordinates(sites.sp))
round(mean(sites.d),2)
round(range(sites.d),2)
#View(sites.on2)
