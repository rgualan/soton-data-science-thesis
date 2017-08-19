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
source("util/my_helper.R")

## Gloval variables ##################################################################
paper <- setupPaper()
numDays <- as.integer((convertStringToPOSIXct("2016-12-31")
                       -convertStringToPOSIXct("2016-01-01"))+1)


## Read data #########################################################################
if(F){
  ## Build the list of active stations starting from the data
  # d <- read.csv("data/epa/epa_daily/2016/daily_44201_2016_ozone.csv", header=T)
  # d$Date <- as.POSIXct(d$Date.Local,format="%Y-%m-%d", tz="GMT")
  # ## Filter California
  # d <- d[d$State.Name=="California", ]
  # d$Site <- as.factor(sprintf("%03d-%04d-%01d",d$County.Code,d$Site.Num,d$POC))
  # d$Station.Code <- as.factor(sprintf("%03d-%04d",d$County.Code,d$Site.Num))
  # ## Generic variable name
  # names(d)[names(d)=="Arithmetic.Mean"] <- "Ozone"
  # ## Keep relevant fields
  # d <- d[,c("Station.Code","Site","POC", "Date","Ozone")]
  # head(d)
  # saveRDS(d, "data/epa/epa_daily/2016/daily_44201_2016_ozone_ca.RDS")
  #d <- readRDS("data/epa/epa_daily/2016/daily_44201_2016_ozone_ca.RDS")
  #sites <- getSites(d)
  #d <- merge(d,sites)  # Implicit removing of stations without Location.Setting
  #View(sites)
  
  
  ## Melt POC ###########################################################################
  ## Stations with several POC
  #table(d[,c("Station.Code","POC")])
  d <- aggregate(Ozone~Station.Code+Date,d,mean)
  
  ## Disable URBAN stations ########################################################
  ## Some papers only use RURAL stations
  ## But, others use the type of station as a covariate
  ## Therefore, we should keep all the stations
  # head(sites)
  # d2 <- merge(d, sites[,c("Station.Code","Location.Setting")])
  # head(d2)
  # levels(d2$Location.Setting)
  # unique(d2$Location.Setting)
  # nrow(d)
  # d <- d2[d2$Location.Setting %in% c("RURAL",""),-4]
  # nrow(d)
  # unique(d$Station.Code)

  ## Disable stations with no enough data #########################################
  recordsByStation <- aggregate(Ozone~Station.Code, d, length)
  names(recordsByStation)[ncol(recordsByStation)] <- "Count"
  #View(recordsByStation);
  minNumRecords <- 0.90*numDays
  ## Disable stations
  ## Before:
  nrow(recordsByStation)
  ## After 
  sites.on <- recordsByStation[recordsByStation$Count>=minNumRecords,]
  nrow(sites.on)
  # table(merge(sites.on,sites_md)$Location.Setting)
  ## Apply filter
  d <- d[d$Station.Code %in% sites.on$Station.Code,]
  sites <- getSites(d)
  nrow(sites)

  ## Save dataset #####################################################################
  # d <- d[order(d$Station.Code,d$Date),]
  # names(d)[3] <- "Ozone" 
  # sites <- sites[order(sites$Station.Code),]
  # saveRDS(d, "data/epa/epa_daily/2016/california_ozone.RDS")
  # saveRDS(sites, "data/epa/epa_daily/2016/california_ozone_sites.RDS")  
}

## Cleaned !!!
d <- readRDS("data/epa/epa_daily/2016/california_ozone_2.RDS")
sites <- getSites(d)

## Initial exploration ##############################################################################
printPlot(paper,"img/eda/hist_ozone.jpeg",5,5,FUN=function(){
  hist(d$Ozone, main="")  
})

if(F){
  ggplot(d) + 
    geom_line(aes(x=Date,y=Ozone,col=Station.Code), alpha=0.5, size=0.5) +
    theme(legend.position="none")
}

## Plot active stations ##################################################################
## Map of California
printPlot(paper,"img/eda/ca_ozone.jpeg",6,6,FUN= function(){
  p <- ggplot(getCAmap()) +
    geom_polygon(aes(x = long, y = lat, group = group), fill = "white", colour = "black") +
    geom_point(data = sites, aes(x = Longitude, y = Latitude, fill=Location.Setting), #, fill=Location.Setting
               alpha = 0.75, shape=24, size=2) +
    labs(x = "Longitude", y = "Latitude", fill="Type") +
    coord_quickmap() +
    theme(legend.justification = c("right", "top"), legend.position = c(.95, .95),
          legend.box.background = element_rect(), legend.box.margin = margin(6, 6, 6, 6))
  plot(p)
})

## Concentration plus Daily mean TS ##################################################################
## Daily particle concentrations for the monitoring sites sorted from the top 
## to the bottom by decreasing longitude (from west to east). 
## The bottom panel shows the daily median concentrations across
## the time-series.
# levelplot(Ozone~Date*Station.Code,
#           dtmp[dtmp$Station.Code %in% sample(unique(dtmp$Station.Code),50),],
#           cuts=10,col.regions=rev(brewer.pal(11,"Spectral")),
#           scales=list(y=list(cex=.6)))
printPlot(paper, "img/eda/heatmap_ca_ozone.jpeg", 6, 7, FUN= function(){
  heatmapPlusTs(d, "Ozone", "Ozone (ppb)")
})

## Check NAs ##########################################################
## Plot a data concentration as matrix (sites x date)
if(F){
  NS<-50
  for(i in 1:ceiling(nrow(sites)/50)){
    print(i)
    print(levelplot(Ozone~Date*Station.Code, 
                    d[d$Station.Code %in% sites$Station.Code[(NS*(i-1)+1):(NS*i)],],
                    cuts=10,col.regions=rev(brewer.pal(11,"Spectral")),
                    scales=list(y=list(cex=.7))))
    if(!par("ask")) readline("Continue?")
  }
}

## Statistics about NAs
recordsByStation <- aggregate(Ozone~Station.Code, d, length)
names(recordsByStation)[2] <- "Count" 
recordsByStation$Missing <- numDays-recordsByStation$Count
recordsByStation$MissingP <- (recordsByStation$Missing/numDays)*100
mean(recordsByStation$MissingP)
range(recordsByStation$MissingP)


## Statistics about Distance
sites.sp <- convertDataToSp(sites)
sites.d <- dist(coordinates(sites.sp))
round(mean(sites.d),2)
round(range(sites.d),2)
#View(sites.on2)


## Compare urban and suburban types ############################################################
## How different are the time series from the different types of location settings
d2 <- merge(d,sites)
aggByType <- aggregate(Ozone~Location.Setting+Date,d2,mean)
#View(aggByType)
if(F){
  ggplot(d2, aes(Date, Ozone, colour=Location.Setting)) + 
    geom_line(size=0.4) +
    labs(x = "Date", y = "Ozone (ppb)") +
    theme(legend.position = "top")
  
  cor(aggByType$Ozone[aggByType$Location.Setting=="RURAL"],  ## These are the most correlated
      aggByType$Ozone[aggByType$Location.Setting=="SUBURBAN"],
      use="pairwise.complete.obs")
  cor(aggByType$Ozone[aggByType$Location.Setting=="URBAN"],
      aggByType$Ozone[aggByType$Location.Setting=="SUBURBAN"],
      use="pairwise.complete.obs")
  cor(aggByType$Ozone[aggByType$Location.Setting=="RURAL"],  ## These are the least correlated
      aggByType$Ozone[aggByType$Location.Setting=="URBAN"],
      use="pairwise.complete.obs")
}

## Individual comparisons 
if(F){
  ## Map plot with labels by station
  ## for comparing how close stations from different types are
  sites$Short.ID <- 1:nrow(sites)
  ggplot(getCAmap()) +
    geom_polygon(aes(x = long, y = lat, group = group), fill = "white", colour = "black") +
    geom_text(data = sites, aes(x = Longitude, y = Latitude, label=Short.ID, 
                                color=Location.Setting),
              alpha = 0.75, size=2, check_overlap = T) +
    labs(x = "Longitude", y = "Latitude", fill="Type") +
    coord_quickmap() +
    theme(legend.position = "top")
  
  #ids <- c(149, 150) # (urban vs suburban)
  #ids <- c(59, 60)
  ids <- c(81, 82)
  sc <- sites$Station.Code[sites$Short.ID %in% ids]
  a <- d[d$Station.Code==sc[1],]
  b <- d[d$Station.Code==sc[2],]
  a <- fillMissingDates(a)
  b <- fillMissingDates(b)
  plot(Ozone~Date,a,type="l")
  lines(Ozone~Date,b,col=2)
  cor(a$Ozone,b$Ozone,use="pairwise.complete.obs")^2
  dist(sites[sites$Station.Code %in% sc,c("Longitude","Latitude")])
  ## NOTES:
  ## 1) 0.9509551 - 0.2200834
  ## 2) 0.4002135 - 0.6695808
  ## 3) 0.8907591 - 0.2822242
  ##  
  ## Sites changed!!!
  ## 1) They are basically the same, even though they belong to different Location types!!!
  ## 2) Weakly correlated 0.4
  ## 3) Strongly correlated
  ##
  ## Distance is a good indicator of how correlated two stations are.
  ## Location.Setting is not
}


## Correlation versus distance #####################
## Create scatterplot showing correlation of daily data for pairs of monitoring sites
## as a function of their separation distance.
if(paper){
  sites <- getSites(d)
  ks <- sites[,c("UTM.X","UTM.Y")] 
  rownames(ks) <- sites[,c("Station.Code")]
  #head(ks)
  m <- as.matrix(dist(ks))
  dim(m)
  m[1:5,1:5]
  ## Build dataframe with distances and correlations
  r2 <- data.frame()
  for(i in 2:nrow(m)){
    for(j in 1:(i-1)){
      s1 <- rownames(m)[i]
      s2 <- colnames(m)[j]
      a <- m[i,j]
      if(a<200){
        ts1 <- d[d$Station.Code==s1,c("Date","Ozone")]
        ts2 <- d[d$Station.Code==s2,c("Date","Ozone")]
        tsC <- merge(ts1,ts2,by="Date")
        b <- cor(tsC[,2],tsC[,3],use="pairwise.complete.obs")^2 
        newRow <- data.frame(S1=s1,S2=s2,Distance=a,Correlation=b)
        r2 <- rbind(r2,newRow)
      }
    }
  }
  ## Plot
  printPlot(paper,"img/eda/correlation_vs_distance_200_ozone.jpeg",5,5,FUN=function(){
    plot(Correlation~Distance,r2, ylab="Correlation", xlab="Distance (Km)")
  })
  printPlot(paper,"img/eda/correlation_vs_distance_50_ozone.jpeg",5,5,FUN=function(){
    plot(Correlation~Distance,r2[r2$Distance<=50,], ylab="Correlation", xlab="Distance (Km)")
  })
}