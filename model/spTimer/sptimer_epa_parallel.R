## demo("nyExample")

## Clean the environment
rm(list=ls())

## Load libraries
library(spTimer)
library(mgcv)
library(spBayes)
library(maps)
library(parallel)
source("util/my_helper.R")

## Global variables ###############################################################
paper <- setupPaper()


## Read data #######################################################################
epa <- readRDS("data/epa/epa_daily/2016/california_ozone_plus_rcov_3.RDS")
epa <- epa[order(epa$Station.Code, epa$Date),]
sites <- getSites(epa)
## Test
epa <- addDoyField(epa)
epa <- addDowField(epa)
## Standardize variable 
epa$sOzone <- scale(epa$Ozone)[,1]

## Build design matrix ###############################################################
## Create additional columns for type factor (method B)
## The Rural factor is the base
## The other types are addicional information ???
epa$fac2 <- 0
epa$SuburbanTemp <- 0
epa$SuburbanRH <- 0
epa$fac3 <- 0
epa$UrbanTemp <- 0
epa$UrbanRH <- 0
epa$fac2[epa$Location.Setting=="SUBURBAN"] <- 1
epa$SuburbanTemp[epa$Location.Setting=="SUBURBAN"] <- epa$Temperature[epa$Location.Setting=="SUBURBAN"]
epa$SuburbanRH[epa$Location.Setting=="SUBURBAN"] <- epa$RH[epa$Location.Setting=="SUBURBAN"]
epa$fac3[epa$Location.Setting=="URBAN"] <- 1
epa$UrbanTemp[epa$Location.Setting=="URBAN"] <- epa$Temperature[epa$Location.Setting=="URBAN"]
epa$UrbanRH[epa$Location.Setting=="URBAN"] <- epa$RH[epa$Location.Setting=="URBAN"]
#aggregate(fac2~Location.Setting,epa,sum)
#epa=epa[,names(epa) != "Location.Setting"]; # Del Location.Setting

covariates <- c("Temperature", "RH", "fac2", "SuburbanTemp", "SuburbanRH", 
                "fac3", "UrbanTemp", "UrbanRH")
round(epa[1:5,covariates],2)


## 10-fold cross validation #################################################################################
#fm <- sOzone ~ 1
fm <- sOzone ~ Temperature+Elevation+Doy
#fm <- sOzone ~ Temperature+RH+sqrtWind+Elevation+Location.Setting+Doy+Dow.name+Dow.number+isWeekday

folds <- readRDS("output/folds.RDS")
cl <- makeCluster(11, outfile="")
clusterExport(cl, c("epa","sites","paper","fm","folds","ticToc"))
tryCatch({
  out <- clusterApply(cl, 1:10, function(k){
    library(spTimer)
    source("util/my_helper.R")
    print(paste("Fold",k))
    
    ## Split data
    epa.train <- epa[epa$Station.Code %in% sites$Station.Code[folds!=k],]
    epa.test <- epa[epa$Station.Code %in% sites$Station.Code[folds==k],]
    epa.test$sOzoneH <- NA  ## For the prediction
    
    ## Fit model
    ticToc({
      post.gp <- spT.Gibbs(formula = fm, data = epa.train, model = "GP", 
                           coords = ~UTM.X+UTM.Y, 
                           spatial.decay = spT.decay(distribution = Gamm(2,1), tuning = 0.1),
                           report=2)
    })
    #print(post.gp)
    summary(post.gp)

    ## Prediction
    ticToc({
      pred.gp <- predict(post.gp, newdata=epa.test, newcoords=~UTM.X+UTM.Y)
    })
    epa.test$sOzoneH <- c(pred.gp$Median)
    ## Fold metrics
    m <- evaluatePredictions(epa.test$sOzone, epa.test$sOzoneH)
    print(sprintf("Fold-%d. Metrics:",k))
    print (m)
    
    return(epa.test)
  })
}, error = function(e){
  print(e)
}, finally = {
  stopCluster(cl)
})

## Save results
saveRDS(out, "output/spTimer/out.p.RDS")
#out <- readRDS("output/spTimer/out.p.RDS")


## Figures 8 (a) -- (d)
## Check convergence
# dev.new()
# par(mfrow = c(2, 2))
# plot(post.gp$betap[1,], type = "l", main = "Intercept", xlab = "Iterations", ylab = "") #, ylim = c(-1, 6)
# plot(post.gp$betap[2, ], type = "l", main = "cMAXTMP", xlab = "Iterations", ylab = "") # , ylim = c(0.06, 0.25)
# # plot(post.gp$betap[3, ], type = "l", main = "WDSP", xlab = "Iterations", ylab = "")
# # plot(post.gp$betap[4, ], type = "l", main = "RH", xlab = "Iterations", ylab = "")
# par(mfrow = c(1, 1))
# readline("Continue?")
# dev.off()


## Calculate metrics ##############################################################
epa.out <- do.call("rbind", out)
## Metrics
metrics <- getMetricsByStationFromDF(epa.out,"sOzone","sOzoneH")
saveRDS(metrics, "output/spTimer/out.metrics.RDS")
