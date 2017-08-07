# Clean environment
rm(list=ls())

# Load libraries
library(raster)
library(ncdf4)


# Function for extracting time series from a nc file ##################################
extractTimeSeriesFromNc <- function(ncFile, dateA, dateB, thePoints){
  # read the netCDF file as a raster layer
  r <- raster(ncFile)
  #r; print(r)
  
  covariate <- data.frame( date=c(),Station.ID=c(),tmp=c())
  #i=1
  for(i in 1:nbands(r)){
    #for(i in 1:5){
    print(i)  
    r <- raster(ncFile, band=i)
    z <- as.POSIXct(getZ(r), format="%Y-%m-%d")
    if(z>=dateA & z<=dateB){
      print("Processing")
      vals <- extract(r, thePoints, 
                      method='bilinear', fun=mean, na.rm=TRUE)
      tmp <- data.frame( Date=rep(z,length(vals)),
                         Station.ID=stations$Station.ID,
                         var=vals)
      covariate <- rbind(covariate, tmp)
    }
  }

  return(covariate)  
}





# Read main df ###############################################################
d <- read.csv("data/FR_AQeReporting_2013-2015/France_data_byday_0.csv")
d$Date <- as.POSIXct(d$Date, format="%Y-%m-%d", tz="UTC")
#View(d)

# Points (stations) to extract data from
stations <- unique(d[,c("Station.ID","Longitude","Latitude")])
dateA <- min(d$Date)
dateB <- max(d$Date)
#dateC <- as.POSIXct("2015-03-01", format="%Y-%m-%d", tz="UTC")
thePoints <- cbind(stations$Longitude,stations$Latitude)

## Extract temperature (TMP) ###################################################
temperature <- extractTimeSeriesFromNc("data/e_obs_grid/tg_0.25deg_reg_v15.0-2015-box.nc", 
                                       dateA, dateB, thePoints)
names(temperature)[which(names(temperature)=="var")] <- "TMP" # TODO: Use short name instead
## Extract RAINFALL (RAIN) ###################################################
rain <- extractTimeSeriesFromNc("data/e_obs_grid/rr_0.25deg_reg_v15.0-2015-box.nc", 
                                       dateA, dateB, thePoints)
names(rain)[which(names(rain)=="var")] <- "RAIN" 
## Extract PRESSURE (SLP) ###################################################
pressure <- extractTimeSeriesFromNc("data/e_obs_grid/pp_0.25deg_reg_v15.0-2015-box.nc", 
                                dateA, dateB, thePoints)
names(pressure)[which(names(pressure)=="var")] <- "PRESS" # TODO: Use short name instead

# Merge dataframes
complete <- d
complete <- merge(complete,temperature,all.y=T)
complete <- merge(complete,rain,all.y=T)
complete <- merge(complete,pressure,all.y=T)

#View(complete)
#dim(complete)

complete <- complete[order(complete$Date,complete$Station.ID),]
complete <- complete[,c("Station.ID","Date","Altitude","Longitude","Latitude","UTMX","UTMY",
                        "TMP","RAIN","PRESS","PM10")]

# Analyze NAs
aggNas <- aggregate(cbind(TMP,RAIN,PRESS,PM10)~Station.ID, FUN= function(x){sum(is.na(x))}, data=complete, na.action=na.pass)
aggNas2 <- aggNas[aggNas$TMP>0 | aggNas$RAIN>0 | aggNas$PRESS>0 | aggNas$PM10>0,]
#View(aggNas2)
# Ignore satations with NAs
complete <- complete[!(complete$Station.ID %in%
                         aggNas$Station.ID[aggNas$TMP>30 | aggNas$RAIN>30 | aggNas$PRESS>30]),]
#dim(complete)

# Write training and validation
#completeA <- complete[complete$Date<dateC,]
#completeB <- complete[complete$Date>=dateC,]
stations_train <- read.csv(file="data/FR_AQeReporting_2013-2015/coordinates.csv", header=T)
stations_val <- read.csv(file="data/FR_AQeReporting_2013-2015/coordinates_val.csv", header=T)
completeA <- complete[complete$Station.ID %in% stations_train$Station.ID,]
completeB <- complete[complete$Station.ID %in% stations_val$Station.ID,]
write.csv(completeA, file="data/FR_AQeReporting_2013-2015/France_data_byday.csv", row.names = F)
write.csv(completeB, file="data/FR_AQeReporting_2013-2015/France_data_byday_val.csv", row.names = F)


