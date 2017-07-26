# Clean environment #################################################################
rm(list=ls())

# Libraries
library(openair)

## Stations ###################################################################
## Read data
load("data/kcl/sites.RData")
nrow(sites) # 900

## Download active stations
## Period
dateA <- as.POSIXct("2015-01-01", format="%Y-%m-%d", tz="UTC")
dateB <- as.POSIXct("2015-04-30", format="%Y-%m-%d", tz="UTC")
## Area
minLon <- -0.567 
maxLon <- 0.312
minLat <- 51.2 
maxLat <- 51.8
## Filter
sites2 <- sites[sites$OpeningDate<=dateA 
                & (sites$ClosingDate>=dateB | is.na(sites$ClosingDate))
                & sites$Longitude >= minLon & sites$Longitude <= maxLon
                & sites$Latitude >= minLat & sites$Latitude <= maxLat,]
nrow(sites2)

# counter <- 1
for(s in sites2$SiteCode[!is.na(sites2$SiteCode)]){
  print(s)
  fileName <- paste0(s,"_2015.RData")
  URL <- paste0("http://www.londonair.org.uk/r_data/",fileName)
  destFile <- paste0("data/kcl/",fileName)
  tryCatch({
    download.file(URL,destFile)
  }, error = function(e){print("Didn't work!")})
  
  ## Limit
  # counter = counter + 1
  # if(counter == 10){break}
}


## Randomly check some of the downloaded files ###################################
files <- list.files("data/kcl/2015")
for(f in files[sample(length(files),5)]){
  print(f)
  load(paste0("data/kcl/2015/",f))
  summaryPlot(x)
  readline("Continue?")
}


