## Clean environment #################################################################
rm(list=ls())

## Libraries
library(sp)
library(lattice)
library(RColorBrewer)
library(ggplot2)
library(maps)
source("util/my_helper.R")

## Settings ##########################################################################
printPlots <- T

## Read data #########################################################################
sites <- read.csv("data/epa/sites/aqs_sites.csv", header=T)


## Plot all ACTIVE  stations (based on Closed.Date field) ############################
sites.active <- sites[sites$Site.Closed.Date=="",]

mapUSA <- getUSAmap()
sites.active <- sites.active[sites.active$Longitude> mapUSA@bbox[1,1]
                             & sites.active$Longitude< mapUSA@bbox[1,2]
                             & sites.active$Latitude> mapUSA@bbox[2,1]
                             & sites.active$Latitude< mapUSA@bbox[2,2],]
sites.active$Location.Setting[sites.active$Location.Setting==""] <- "UNKNOWN"

if(printPlots) jpeg("img/eda/all_sites.jpeg", 8, 4, "in", bg="white", res=150)
ggplot(mapUSA) +
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", colour = "black") +
  geom_point(data = sites.active, 
             aes(x = Longitude, y = Latitude, colour=Location.Setting), 
             alpha = .75, size=0.75, shape=17) + #
  labs(x = "Longitude", y = "Latitude", shape="Type") +
  coord_quickmap() + 
  theme(legend.position = "top")
if(printPlots) dev.off()


## Select California ###########################################################################
sites <- sites[sites$State.Name=="California",]


## Set Station.Code and relevant columns ###########################################################################
sites$Station.Code <- as.factor(sprintf("%03d-%04d",sites$County.Code,sites$Site.Num))
sites <- sites[,c("Station.Code","State.Code", "County.Code", "Site.Number", 
                        "Latitude", "Longitude", "Datum", "Elevation", "Location.Setting",
                        "State.Name")] # Keep relevant fields only
#View(sites)


## Obtain locations in kilometers using UTM  ####################################################
sites.sp <- sites
coordinates(sites.sp) <- ~Longitude+Latitude
proj4string(sites.sp) <- "+proj=longlat +datum=WGS84"
utm11 = CRS("+proj=utm +zone=11 +datum=WGS84 +ellps=WGS84 +units=km")
sites.sp <- spTransform(sites.sp,utm11) 
sites_utm <- as.data.frame(sites.sp)
names(sites_utm)[9:10] <- c("UTM.X", "UTM.Y")
sites <- merge(sites,sites_utm[,c("Station.Code","UTM.X", "UTM.Y")])
#View(sites)


## Rename URBAN level ###############################################################
levels(sites$Location.Setting)
aggregate(Station.Code~Location.Setting,sites,length)
levels(sites$Location.Setting)[5] <- "URBAN"  ## Make the name shorter



# ## Combine Suburban and Urban sites (Bayesian2014) ###########################
# ## Combine URBAN and Suburban -> Urban
# sites$Location.Setting2 <- sites$Location.Setting
# sites$Location.Setting2[sites$Location.Setting2=="SUBURBAN"]<-"URBAN"
# aggregate(Station.Code~Location.Setting2,sites,length)
# ## What are the unknown sites
# mapCA <- mapUSA[mapUSA$NAME_1=="California",]; proj4string(mapCA)
# ggplot(mapCA) +
#   geom_polygon(aes(x = long, y = lat, group = group), fill = "white", colour = "black") +
#   geom_point(data = sites[sites$Location.Setting2=="UNKNOWN",], 
#              aes(x = Longitude, y = Latitude, colour=Location.Setting), 
#              alpha = .75, size=0.75, shape=17) + #
#   labs(x = "Longitude", y = "Latitude", shape="Type") +
#   coord_quickmap() + 
#   theme(legend.position = "top")
# ## Assign the UNKNOWN sites to the Type of the closest station
# sites[sites$Location.Setting2=="UNKNOWN",]
# #View(sites[sites$Location.Setting=="UNKNOWN",])
# library(FNN)
# ks <- sites[! sites$Location.Setting2 %in% c("UNKNOWN",""),c("Station.Code","UTM.X","UTM.Y")] # known sites
# 
# for(s in sites$Station.Code[sites$Location.Setting %in% c("UNKNOWN","")]){
#   print(s)
#   ## Get closest neighbour
#   print(knnx.dist(ks[,-1], sites[sites$Station.Code==s,c("UTM.X","UTM.Y")], k=1))
#   i <- knnx.index(ks[,-1], sites[sites$Station.Code==s,c("UTM.X","UTM.Y")], k=1)
#   s_ <- ks[i,1]
#   sites$Location.Setting2[sites$Station.Code==s] <-
#     sites$Location.Setting2[sites$Station.Code==s_]
# 
#   if(F){
#     p<-ggplot(sites) +
#       geom_polygon(aes(x = long, y = lat, group = group), data=mapCA, fill = "white", colour = "black") +
#       geom_point(aes(x = Longitude, y = Latitude, colour=Location.Setting2),
#                  alpha = .50, size=0.75, shape=17) +
#       geom_point(data = sites[sites$Station.Code==s,],
#                  aes(x = Longitude, y = Latitude), 
#                  alpha = .75, size=2, shape=24, fill="red", colour="black") + #
#       geom_point(data = sites[sites$Station.Code==s_,], 
#                  aes(x = Longitude, y = Latitude, fill=Location.Setting2), 
#                  alpha = .75, size=2, shape=24, colour="black") + #
#       labs(x = "Longitude", y = "Latitude", shape="Type") +
#       coord_quickmap() + 
#       theme(legend.position = "top")
#     plot(p)
#     readline("Continue?")
#   }
# }
# ## Before
# aggregate(Station.Code~Location.Setting,sites,length)
# ## After
# aggregate(Station.Code~Location.Setting2,sites,length)
# sites$Location.Setting <- droplevels(sites$Location.Setting2)
# #plot(sites$Location.Setting)
# sites<-sites[,-length(names(sites))]
# #head(sites)


## How to fill the missing Location.Setting in some stations??? ################################
## Define the type of a station based on the type of the closest neighbour 
## This is a naive impolementation
## There are better options
if(F){
  ## Stations without Location.Setting
  ggplot(getCAmap()) +
    geom_polygon(aes(x = long, y = lat, group = group), fill = "white", colour = "black") +
    geom_point(data = sites, 
               aes(x = Longitude, y = Latitude, colour=Location.Setting), 
               alpha = .75, size=0.75, shape=17) + 
    geom_point(data = sites[sites$Location.Setting
                            %in% c("UNKNOWN",""),], 
               aes(x = Longitude, y = Latitude), 
               alpha = .75, size=2, shape=17, colour="red") + #
    labs(x = "Longitude", y = "Latitude", shape="Type") +
    coord_quickmap() + 
    theme(legend.position = "top")
}
## Notes:
## This approach could lead to bias in the time series by Location.Setting
## Since the amount of available stations is considerable, 
## the simpler method is to ignore these stations
print("Ignoring stations with no Locattion.Setting information")
nrow(sites[sites$Location.Setting %in% c("UNKNOWN",""),])
sites <- sites[!sites$Location.Setting %in% c("UNKNOWN",""),]




## Save dataset #####################################################################
sites <- sites[order(sites$Station.Code),]
saveRDS(sites, file="data/epa/sites/aqs_sites.RDS")
