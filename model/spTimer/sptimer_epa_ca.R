## demo("nyExample")

## Clean the environment
rm(list=ls())

## Load libraries
library(spTimer)
#library(akima)
#library(coda)
#library(spacetime)
#library(fields)
#library(forecast)
#library(MASS)
library(mgcv)
library(spBayes)
#library(colorspace) 
library(maps)
#library(MBA)
#library(openair)
source("util/my_helper.R")


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


## Split data ####################################################################
folds <- readRDS("output/folds.RDS")
## Fold(1)
epa.train <- epa[epa$Station.Code %in% sites$Station.Code[folds!=1],] 
epa.test <- epa[epa$Station.Code %in% sites$Station.Code[folds==1],]
epa.test$sOzoneH <- NA
#View(epa.train); View(epa.test)



## Simple LM (test) ##############################################################
set.seed(11)
f1 <- sOzone~Temperature+RH+UTM.X*UTM.Y+j+wd+Location.Setting # 0.4214 
f2 <- sOzone~Temperature+UTM.X*UTM.Y # 0.3622
f3 <- sOzone~RH+UTM.X*UTM.Y # 0.2978
lm.fit <- lm(f3, data=epa.train)
summary(lm.fit)
lm.pred <- predict(lm.fit, epa.test, interval="prediction")
evaluatePredictions(epa.test$sOzone, lm.pred[,1])
plot(epa.test$sOzone,type="l"); lines(lm.pred[,1],col=2)

## Simple GAM model (test) ########################################################
gam.fit <- gam(sOzone ~ s(RH), # + s(RH) + s(UTM.X, UTM.Y) + s(cMAXTMP) + s(WDSP) + s(RH) +  # , k = 10
               data = epa.train)
gam.pred <- predict(gam.fit, epa.test, interval="prediction")
evaluatePredictions(epa.test$sOzone, gam.pred)
plot(epa.test$sOzone,type="l"); lines(gam.pred,col=2)
## Notes
## This simple models were used to assess which of the predictors was the most relevant
## Temperature seems to be the most relevant predictor


## Pre-process data ######################################################################


## Fit GP model ######################################################################
set.seed(11)
ticToc({
  post.gp <- spT.Gibbs(formula = sOzone ~ 1, data = epa.train, model = "GP", 
                       coords = ~UTM.X+UTM.Y, #scale.transform = "SQRT", 
                       spatial.decay = spT.decay(distribution = Gamm(2,1), tuning = 0.1),
                       report=10)
})
#saveRDS(post.gp, file="model/spTimer/post.gp.RDS")
#post.gp<-readRDS("model/spTimer/post.gp.RDS")
print(post.gp)
summary(post.gp)
# plot(post.gp) #It works with dev.new()

## Spatial prediction for the GP model
pred.gp <- predict(post.gp, newdata=epa.test, newcoords=~UTM.X+UTM.Y)
print(pred.gp)
names(pred.gp)

## Validation criteria
evaluatePredictions(epa.test$sOzone, c(pred.gp$Median))


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



