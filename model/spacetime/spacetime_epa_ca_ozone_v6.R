## Ref: 
## Spatio-Temporal data in R
## https://www.r-bloggers.com/spatio-temporal-kriging-in-r/
## https://cran.r-project.org/web/packages/gstat/vignettes/st.pdf

## Clean environment
rm(list=ls())

## Libraries
library(gstat)
library(spacetime)
#library(rgdal) ## Not available in the cluster!
library(reshape2)
library(parallel)
source("util/my_helper.R")
source("model/spacetime/testStKriging.R")

## Global variables
set.seed(123)
smplDays <- sort(sample(365,8))
forceRun <- T
paper <- setupPaper()
debugLevel <- T

## Read data #################################################################
epa <- readRDS("data/epa/epa_daily/2016/california_ozone_plus_rcov_3.RDS")
## Add date covariates
epa <- addDoyField(epa)
epa <- addDowField(epa)
## Sites
sites <- getSites(epa)
## Scale target variable #####################################################
epa$sOzone <- scale(epa$Ozone)[,1]  ## scale returns a matrix

## Assemble STFDF ############################################################
## Space dimension
epa.sp <- getSites(epa)
epa.sp <- epa.sp[,c("Station.Code","Latitude","Longitude","Location.Setting","UTM.X","UTM.Y")]
rownames(epa.sp) <- epa.sp$Station.Code
head(epa.sp)
coordinates(epa.sp) <- ~UTM.X+UTM.Y
proj4string(epa.sp) <- getUTMproj()

## Time dimension
epa.tm <- sort(unique(epa$Date))
epa.tm <- as.Date(epa.tm)  ## Ignore time data?? >> Corrects x labels problem in acf
# Combine the objects spatial, data-frame and time-dim into a STIDF:
epa.st <- STFDF(epa.sp,epa.tm,epa[,c("Date","Ozone","sOzone","Temperature",
                                     "RH","Rain","logRain","Wind","sqrtWind","Elevation")]) 
# summary(epa.st)
# stplot(epa.st[,"2016-01-01::2016-01-08","sOzone"])
# dim(epa.st)

## Run simulations in parallel ###############################################################
cl <- makeCluster(4, outfile="")
timeLags <- list(a=0:5, b=0:7, c=0:10)
clusterExport(cl, c("epa.st","paper","timeLags"))
tryCatch({
  clusterApply(cl, 1:length(timeLags), function(tlIndex){
    source("model/spacetime/testStKriging.R")
    print(tlIndex)
    fm1 <- sOzone~1
    fm2 <- sOzone~Temperature
    #out <- test_st_kriging_using_10cv(fm, epa.st, timeLags[[tlIndex]], paper, as.character(tlIndex))
    try({
      out <- test_st_kriging_using_1fold(fm1, epa.st, timeLags[[tlIndex]], paper, as.character(tlIndex))
      saveRDS(out, paste0("output/spacetime/out1_",tlIndex,".RDS"))
    })
    try({
      out <- test_st_kriging_using_1fold(fm2, epa.st, timeLags[[tlIndex]], paper, as.character(tlIndex))
      saveRDS(out, paste0("output/spacetime/out2_",tlIndex,".RDS"))
    })
  })
}, error = function(e){
  print(e)
}, finally = {
  stopCluster(cl)
})
