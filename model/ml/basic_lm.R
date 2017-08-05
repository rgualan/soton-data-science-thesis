## Clean environment
rm(list=ls())

## Load libraries
library(spTimer)

## Global parameters ################################################
#.pardefault
applyGLFilter<-F


## Reading the data files #######################################
spcheck1<-read.table("data/aurn/AURN_data_07_11.txt",header=T)

## Create a dataframe with the following columns:
## index, lon, lat, year, month, day, type, sqrtaqm, obs 
spcheck = spcheck1[,c("index","lon","lat","year","month","day")]
spcheck$date <- as.POSIXct(
  sprintf("%04d-%02d-%02d", spcheck1$year,spcheck1$month,spcheck1$day), tz="GMT")
spcheck$type = spcheck1$type
spcheck$sqrtaqm = sqrt(spcheck1$aqum_no2) ## SQRT transform
spcheck$obs=spcheck1$obs_no2 # the variable! 
spcheck$sqrtobs=sqrt(spcheck1$obs_no2) # the variable! 
## Order by index, date (year, month, day)
spcheck = spcheck[order(spcheck$index, spcheck$date),]
## Check
head(spcheck)


######## Choosing the fit and validation sites ########
sites <- unique(spcheck[,c("index","lon","lat")]) #144

set.seed(11)
sites.val <- sample(sites$index, 44) 
sites.fit <- sites$index[!sites$index %in% sites.val]

## Quick scatter plot 
plot(lat~lon,sites[sites$index %in% sites.fit, ], col=1)
points(lat~lon,sites[sites$index %in% sites.val, ], col=2)

## Validation dataset
prediction_file <- merge(data.frame(index=sites.val), spcheck)
prediction_file <- prediction_file[order(prediction_file$index,prediction_file$date),]
## Training dataset 
fitting_file<-merge(data.frame(index=sites.fit), spcheck)
fitting_file<-fitting_file[order(fitting_file$index,fitting_file$date),]



## Simple LM ######################################################################
lm.1 <- lm(sqrtobs~sqrtaqm+type, fitting_file)
# par(mfrow=c(2,2)); plot(lm.1); par(.pardefault)
# summary(lm.1)
# 

## Simple plot
# d <- spcheck[spcheck$index==4,]
# plot(sqrtobs~date, d, type="l")
# d$yh <- predict(lm.1,newdata = d)
# lines(yh~date, d, col=2)


## Validation
prediction_file$lm.1 <- predict(lm.1, newdata = prediction_file)^2
spT.validation(prediction_file$obs, prediction_file$lm.1)
cor(prediction_file$obs, prediction_file$lm.1, use = "pairwise.complete.obs")^2




