## Clean environment #################################################################
rm(list=ls())

## Libraries
library(sp)
library(lattice)
library(RColorBrewer)
library(ggplot2)
source("util/my_helper.R")

## Settings ##########################################################################
printPlots <- T

## Read data #########################################################################
## Build the list of active stations starting from the data
## Read data
# d <- read.csv("data/epa/epa_daily/2016/daily_TEMP_2016.csv", header=T)
# #View(d)
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
# saveRDS(d, "data/epa/epa_daily/2016/daily_TEMP_2016_ca.RDS")
d <- readRDS("data/epa/epa_daily/2016/daily_TEMP_2016_ca.RDS")

## Relevant period
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


## Melt POC ###########################################################################
## Stations with several POC
table(d[,c("Station.Code","POC")])

dim(d)
d <- aggregate(Measurement~Station.Code+Date,d,mean)
dim(d)


## Initial exploration ##############################################################################
hist(d$Measurement)
ggplot(d) + 
  geom_line(aes(x=Date,y=Measurement,col=Station.Code), alpha=0.5, size=0.5) +
  theme(legend.position="none")



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
  #readline("Continue?")
}



## Disable stations with no enough data #########################################
recordsByStation <- aggregate(Measurement~Station.Code, d, length)
names(recordsByStation)[ncol(recordsByStation)] <- "Count"
#View(recordsByStation);
minNumRecords <- 0.80*as.integer((dateB-dateA)+1)
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

if(printPlots) jpeg("img/eda/ca_temperature.jpeg", 6, 6, "in", bg="white", res=150)
ggplot(mapCA) +
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", colour = "black") +
  geom_point(data = sitesDs, aes(x = Longitude, y = Latitude, fill=Location.Setting),
             alpha = .75, shape=21, size=2) +
  labs(x = "Longitude", y = "Latitude", fill="Type") +
  coord_quickmap() + 
  theme(legend.justification = c("right", "top"), legend.position = c(.95, .95), 
        legend.box.background = element_rect(), legend.box.margin = margin(6, 6, 6, 6))
if(printPlots) dev.off()


## Check NAs and Time series ##################################################################
## Concentration map
print(levelplot(Measurement~Date*Station.Code,d,
                cuts=10,col.regions=rev(brewer.pal(11,"Spectral")),
                scales=list(y=list(cex=.7))))

heatmapPlusTs(d, "Temperature (°F)", "img/eda/heatmap_ca_temperature.jpeg")




## Notes:
## Maybe remove the station with the event! 
range(d$Measurement)
hist(d$Measurement)
max(d$Measurement)


## Save dataset #####################################################################
d <- d[order(d$Station.Code,d$Date),]
names(d)[3] <- "Temperature" 
saveRDS(d, "data/epa/epa_daily/2016/california_temperature.RDS")
saveRDS(sitesDs, "data/epa/epa_daily/2016/california_temperature_sites.RDS")

