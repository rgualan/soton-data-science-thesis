## Clean environment
rm(list=ls())

## Load required packages ###
library(fields)
library(raster)
library(spatial.tools)
library(gdalUtils)
library(rgdal)
library(gstat)
library(automap)
library(imputation)
library(spTimer)
library(kknn)
source("util/my_helper.R")

## Global variables
paper <- setupPaper()

## Read data ##############################################################################
epa <- readRDS("data/epa/epa_daily/2016/california_ozone_plus_cov.RDS")
## Add date covariates
epa <- addDoyField(epa)
epa <- addDowField(epa)
## Sites
sites <- getSites(epa)
#str(epa)

## Scale target variable
epa$sOzone <- scale(epa$Ozone)[,1]  ## scale returns a matrix

## Use stations with data in the chosen covariates
covByStation <- aggregate(cbind(Ozone,Temperature)~Station.Code,epa,length)
covByStation <- covByStation[covByStation$Ozone>350 & covByStation$Temperature>350,]
epa <- epa[epa$Station.Code %in% covByStation$Station.Code,]
#apply(epa,2,function(x){sum(is.na(x))})
## NOTE
## WARNING: Count the amount of data of the variables 
## Wind speed and RH, in case they are used as covariates

## Plot "semi-complete" stations and test station

## Chose a test station for injecting missing data
## For the slightly isolated by central position: 107-0009 
# sites <- getSites(epa2)
# sites[sites$Longitude> -119 & sites$Longitude< -118.5 & sites$Latitude> 36 & sites$Latitude<37.5, ]
testStation <- "107-0009"

plotStations(paper, covByStation$Station.Code, 
             "img/preprocessing/idw/complete_stations.jpeg",6,6, testStation)

## Inject a random block of missing data and sparse missing data ############### 
## and some random missing data
## sOzone : original (intact)
## sOzone2 : original with injected missing data (predicted)
epa$sOzone2 <- epa$sOzone
set.seed(1345)
randomDay <- sample(getStudyDays(), 1)
epa$sOzone2[epa$Station.Code == testStation &
               epa$Date>=randomDay & epa$Date<=(randomDay+15*24*60*60)] <- NA
rds <- sample(getStudyDays(), 10)
for(i in 1:length(rds)){
  epa$sOzone2[epa$Station.Code == testStation &
                 epa$Date>=rds[i] & epa$Date<=(rds[i]+2*24*60*60)] <- NA
}
#plot(sOzone2~Date,epa[epa$Station.Code==testStation,],type="l")


## GBM ####################################################################################
## Design Matrix
#apply(epa[epa$Station.Code==testStation,],2,function(x){sum(is.na(x))})
#apply(epa[epa$Station.Code!=testStation,],2,function(x){sum(is.na(x))})
epaGbm <- epa[!is.na(epa$Temperature),]
X <- epaGbm[,c("sOzone2","Temperature","UTM.X","UTM.Y","Elevation","Location.Setting","Doy","Dow")]
#apply(X,2,function(x){sum(is.na(x))})  

## Generalized Boosted Regression 
ticToc(
  gbm.fit <- gbmImpute(X, max.iters = 2, cv.fold = 5, verbose=T)
)
epaGbm$sOzone3 <- gbm.fit$x$sOzone2[,1]  ## The output is a matrix
# dim(epaGbm[epaGbm$Station.Code==testStation,])
# plot(sOzone~Date,epaGbm[epaGbm$Station.Code==testStation,],type="l")
# lines(sOzone3~Date,epaGbm[epaGbm$Station.Code==testStation,],col=2)

## Gaussian Process ######################################################################
## First, it is necessary to fill the missing values in Temperature for training
epaGp <- epa
## Simple test
if(F){
  epaGp$Temperature[epaGp$Station.Code==testStation & 
                     epaGp$Date>="2016-05-01" & epaGp$Date<="2016-05-15"] <- NA
  a <- imputeTS::na.kalman(epaGp$Temperature[epaGp$Station.Code==testStation])  ## Simple imputation
  a <- imputeTS::na.ma(epaGp$Temperature[epaGp$Station.Code==testStation])  ## Simple imputation
  plot(epaGp[epaGp$Station.Code==testStation,]$Date,a,col=2,type="l")
  lines(Temperature~Date,epaGp[epaGp$Station.Code==testStation,],type="l",lwd=2)
}
epaGp$Temperature <- imputeTS::na.ma(epaGp$Temperature)

ticToc(
  simpleGp <- spT.Gibbs(
    #formula = sOzone2~Temperature+Elevation+Location.Setting+Doy+Dow,
    formula = sOzone2~Temperature+Elevation+Doy+Dow,
    model = "GP",
    data = epaGp[epaGp$Station.Code!=testStation,], 
    coords = ~Longitude + Latitude, #scale.transform = "SQRT",
    #newdata = epaGp[is.na(epaGp$sOzone2),],
    #newcoords = ~Longitude + Latitude,
    #time.data = spT.time(366),
    spatial.decay = spT.decay(distribution = Gamm(2, 1), tuning = 0.1))
)
simpleGp.pred <- predict(simpleGp, newdata=epaGp[epaGp$Station.Code==testStation,], 
                   newcoords = ~Longitude + Latitude)

## Kernel K-nearest neighbours ##################################################################
# Option a) incomplete
# When using an algorithm where the output depends on distance calculation (as is the case in 
# k-nearest-neighbors) it is recommended to first scale the data
# Option b) functional
# This nearest neighbor method expands knn in several directions.  First it can be used not only for
# classification,  but also for regression and ordinal classification.   Second it uses kernel functions
# to weight the neighbors according to their distances.  In fact, not only kernel functions but every
# monotonic decreasing will work fine.
# Is used in rattle
epaKnn <- epa[,c("Station.Code","sOzone2","Temperature","Elevation","Doy","Dow")]
#epaKnn[,-1] <- scale(epaKnn[,-1])
fit.kknn <- kknn(sOzone2 ~ Temperature+Elevation+Doy+Dow, 
                 epaKnn[epaKnn$Station.Code!=testStation & !is.na(epaKnn$sOzone),], 
                 epaKnn[epaKnn$Station.Code==testStation,])
epaKnn$sOzone5[epaKnn$Station.Code==testStation] <- fit.kknn$fitted.values


# Generated with rattle
# Analyze
# randomForest(formula = sOzone2 ~ .,
#              data = crs$dataset[, c(crs$input, crs$target)],
#              ntree = 2000, mtry = 3, importance = TRUE, replace = FALSE, na.action = randomForest::na.roughfix)


## Test other methods
# fit2 <- SVDImpute(X2, k = 10, num.iters = 10, verbose=F)
# fit2 <- SVTImpute(X2, lambda = 0.1, verbose=F)
# fit2 <- kNNImpute(X, k = 3, verbose = F)
## These methods did not work
## The librayr is not being mantained any more

## Check that the missing data was replaced
# sum(is.na(fit$x)); head(epaGp); head(fit$x)

## Assemble dataframe again ############################################################
test <- epaGbm[epaGbm$Station.Code==testStation,]
test$sOzone4 <- simpleGp.pred$Mean
test$sOzone5 <- epaKnn$sOzone5[epaKnn$Station.Code==testStation]
test2 <- rbind(data.frame(Date=test$Date,Ozone=test$sOzone,Type="Original", Flag=is.na(test$sOzone2)),
               data.frame(Date=test$Date,Ozone=test$sOzone3,Type="GBM", Flag=F),
               data.frame(Date=test$Date,Ozone=test$sOzone4,Type="GP", Flag=F),
               data.frame(Date=test$Date,Ozone=test$sOzone5,Type="KKNN", Flag=F))

printPlot(paper, "img/preprocessing/idw/ts_ozone_others.jpeg",7,3, FUN=function(){ 
  p<-ggplot(test2, aes(x=Date, y=Ozone, colour=Type)) + 
    annotate("rect",
             xmin=test2$Date[test2$Flag]-1*24*60*60,
             xmax=test2$Date[test2$Flag]+1*24*60*60,
             ymin=-Inf, ymax=Inf, alpha=0.75, fill="lightyellow") +
    geom_line() + 
    theme(legend.justification = c("top")) + 
    labs(y="Scaled(Ozone)")
  print(p)
})

# Calculate 
mGbm <- evaluatePredictions(test$sOzone[is.na(test$sOzone2)], test$sOzone3[is.na(test$sOzone2)])
mGp <- evaluatePredictions(test$sOzone[is.na(test$sOzone2)], test$sOzone4[is.na(test$sOzone2)])
mKknn <- evaluatePredictions(test$sOzone[is.na(test$sOzone2)], test$sOzone5[is.na(test$sOzone2)])
metrics <- rbind(mGbm,mGp,mKknn)
rownames(metrics) <- c("gbm","gp","kknn")
metrics

## Notes:
## Imputation is not Interpolation
## So, it shall be tested in a different way 
## Not trying to model a whole time series of an station, but only missing portions 

