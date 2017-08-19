## Clean environment #################################################################
rm(list=ls())

## Libraries
library(sp)
library(lattice)
library(RColorBrewer)
library(ggplot2)
library(grid)
library(gridExtra)
library(dplyr)
library(FNN)
source("util/my_helper.R")

## Read data #########################################################################
d <- readRDS("data/epa/epa_daily/2016/california_ozone.RDS")
smd <- readRDS("data/epa/sites/aqs_sites.RDS")
sites <- unique(d$Station.Code) 
sites <- merge(data.frame(Station.Code=sites),smd)
head(d)
head(sites)
dim(sites)


## Visual help
plotBadStations <- function(sites, theBadStations){
  ggplot() + 
    geom_polygon(aes(x = long, y = lat, group = group), data=getCAmap(), fill = "white", colour = "black") +
    geom_point(data=sites, aes(x = Longitude, y = Latitude, colour=Location.Setting),
               alpha = .50, size=0.75, shape=17) +
    geom_point(data = sites[sites$Station.Code %in% theBadStations,],
               aes(x = Longitude, y = Latitude), 
               alpha = .75, size=2, shape=24, fill="red", colour="black") + #
    labs(x = "Longitude", y = "Latitude", shape="Type") +
    coord_quickmap() + 
    theme(legend.position = "top")
}


## Check close stations and assess correlation 
detectBadStations<-function(d, maxDistance, threshold, minCount){
  sites <- getSites(d)
  ks0 <- sites[,c("Station.Code","UTM.X","UTM.Y")] # known sites
  ks <- ks0[,-1]
  rownames(ks) <- ks0$Station.Code
  #head(ks)
  m <- as.matrix(dist(ks))
  #dim(m)
  #m[1:5,1:5]

  badStations <- c()
  badStationsR2s <- list()
  for(s in unique(d$Station.Code)){
    #s <- sites$Station.Code[1]
    print(paste("Processing",s))
    tsA <- d[d$Station.Code==s,]
    tsA <- fillMissingDates(tsA)
    #View(tsA)
    
    ## Get nearest neighbour
    actualRow <- m[s,]
    kIds <- names(actualRow)[actualRow>0 & actualRow<maxDistance]
    kIds <- kIds[ kIds %in% unique(d$Station.Code) ]
    count <- 0
    rs <- c()
    for(s2 in kIds){
      #s2 <- kIds[1]
      tsB <- d[d$Station.Code==s2,]
      tsB <- fillMissingDates(tsB)
      #View(tsB)
      r <- cor(tsA$Ozone,tsB$Ozone,use="pairwise.complete.obs")
      if (r < threshold){
        count <- count + 1
        rs<-c(rs,r)
        # print(paste("With:",s2,"R:",r))
        # plot(tsA$Ozone,type="l",ylim=range(c(tsA$Ozone,tsB$Ozone)))
        # lines(tsB$Ozone, col=2)
        # readline("Continue?")
      }
    }  
    if(count>minCount){
      print(paste("Incidents: ", count, "R2:", paste(rs,collapse = ",")))
      # if(sites$Location.Setting[sites$Station.Code==s]=="URBAN")
      badStations <- c(badStations,s)
      badStationsR2s[[s]]<-rs
      # else
      #   badStationsR <- c(badStationsR,s)
    }
  }
  output<-list(Stations=badStations, R2s=badStationsR2s)
  return(output)
}

## Approach 1 
## All versus all
if(F){
  badStations1<-detectBadStations(d,50,0.55,2)  
  badStations2<-detectBadStations(d,100,0.4,3)
  theBadStations<-c(badStations1$Stations, badStations2$Stations)
  plotBadStations(sites, theBadStations)
  
  ## Manually check/remove:
  sites[sites$Station.Code %in% theBadStations,c("Station.Code","Elevation","Location.Setting")]
  #theBadStations <- c(theBadStations, "069-0003", "001-2005")
  #theBadStations <- c(theBadStations, badStations)
  #theBadStations <- c(theBadStations, badStationsR)
  d <- d[!d$Station.Code %in% theBadStations,]
  length(unique(d$Station.Code))
  # Save it
  theBadStations <- unique(theBadStations)
  length(theBadStations)
  #saveRDS(theBadStations, file="data/tmp/theBadStations.RDS")
  ## Recover theBadStations
  #theBadStations<-readRDS("data/tmp/theBadStations.RDS")
  saveRDS(d, file="data/epa/epa_daily/2016/california_ozone_2.RDS")
  
  # ## Recover theBadStations
  # theBadStations<-readRDS("data/tmp/theBadStations.RDS")
  # d<-d[!d$Station.Code %in% theBadStations,]
  # saveRDS(d, file="data/epa/epa_daily/2016/california_ozone_2.RDS")
  # ss <- sort(unique(d$Station.Code))
  # sites <- merge(data.frame(Station.Code=ss), smd)
  # # head(sites)
  # # dim(sites)
  # saveRDS(sites, file="data/epa/epa_daily/2016/california_ozone_sites_2.RDS")
}

## Assess Pearson correlation metric between  
## "suspicious" station ts and 
## the mean ts of k closest neighbours to that station
detectBadStationsByTrend<-function(d, maxDistance, threshold, plotDetections=F){
  ## Calculate distance matrix
  sites <- getSites(d)
  ks0 <- sites[,c("Station.Code","UTM.X","UTM.Y")] # known sites
  ks <- ks0[,-1]
  rownames(ks) <- ks0$Station.Code
  #head(ks)
  m <- as.matrix(dist(ks))
  #dim(m)
  #m[1:5,1:5]

  ## Iterate over the stations  
  output <- data.frame()
  for(s in unique(d$Station.Code)){
    ## Current station's ts
    print(paste("Processing",s))
    tsA <- d[d$Station.Code==s,]
    tsA <- fillMissingDates(tsA)
    #View(tsA)
    
    ## Get nearest neighbourS
    actualRow <- m[s,]
    kIds <- names(actualRow)[actualRow>0 & actualRow<maxDistance]
    kIds <- kIds[ kIds %in% unique(d$Station.Code) ]

    ## Calculate the trend
    if(length(kIds)>0){
      trend <- aggregate(Ozone~Date,d[d$Station.Code %in% kIds,],mean)
      trend <- fillMissingDates(trend)
      ## Compare and save
      r <- cor(tsA$Ozone,trend$Ozone,use="pairwise.complete.obs")^2
      output <- rbind(output, data.frame(Station=s, r2=r))
      if(r<threshold & plotDetections){
        print(r)
        plot(Ozone~Date, d[d$Station.Code %in% c(s,kIds),], col=0)
        for(o in kIds){
          lines(Ozone~Date, d[d$Station.Code==o,], col="gray")
        }
        lines(Ozone~Date, tsA, col="red", type="l",lwd=2)
        readline("Continue?")
      }
    }else{
      output <- rbind(output, data.frame(Station=s, r2=NA))
    }
    
    #if(s=="083-4003") stop("Debug!")
  }
  
  ## Output contains NAs
  ## Ignore them before return the bad stations
  output <- output[!is.na(output$r2),]
  
  return(output[output$r2<threshold,])
}




## Approach 2
## One versus trend
if(F){
  badStations1<-detectBadStationsByTrend(d,50, 0.40,F)
  badStations2<-detectBadStationsByTrend(d,100,0.37,F)
  #badStations3<-detectBadStationsByTrend(d,200,0.10,F)
  theBadStationsDF<-rbind(badStations1, badStations2)
  theBadStations <- theBadStationsDF$Station
  plotBadStations(sites, theBadStations)
  
  ## Manually check/remove:
  sites[sites$Station.Code %in% theBadStations,c("Station.Code","Elevation","Location.Setting")]
  #theBadStations <- c(theBadStations, "069-0003", "001-2005")
  #theBadStations <- c(theBadStations, badStations)
  #theBadStations <- c(theBadStations, badStationsR)
  d <- d[!d$Station.Code %in% theBadStations,]
  length(unique(d$Station.Code))
  # Save it
  theBadStations <- unique(theBadStations)
  length(theBadStations)
  #saveRDS(theBadStations, file="data/tmp/theBadStations.RDS")
  ## Recover theBadStations
  #theBadStations<-readRDS("data/tmp/theBadStations.RDS")
  saveRDS(d, file="data/epa/epa_daily/2016/california_ozone_2.RDS")
}


## Check active stations ###########################
ss <- unique(d$Station.Code)
ozone.sites <- merge(data.frame(Station.Code=ss), sites)
head(ozone.sites)
ggplot(sites) +
  geom_polygon(aes(x = long, y = lat, group = group), data=getCAmap(), fill = "white", colour = "black") +
  geom_point(aes(x = Longitude, y = Latitude, colour=Location.Setting),
             alpha = .50, size=0.75, shape=17) +
  labs(x = "Longitude", y = "Latitude", shape="Type") +
  coord_quickmap() + 
  theme(legend.position = "top")
#View(d)

dataByStation <- aggregate(Ozone~Station.Code,d,length)
#View(dataByStation)
dataByStation$pData <- dataByStation$Ozone/366
