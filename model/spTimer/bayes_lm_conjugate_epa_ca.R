## Paper [spTimer]
## ?bayesLMConjugate

## Clean the environment
rm(list=ls())

## Load libraries
library(spTimer)
library(akima)
library(coda)
library(spacetime)
library(fields)
#library(forecast)
library(MASS)
library(mgcv)
library(spBayes)
library(colorspace) 
library(maps)
library(MBA)
source("util/my_helper.R")

## Global parameters ###############################################################
paper <- setupPaper()
forceRun <- T

## Read data #######################################################################
epa <- readRDS("data/epa/epa_daily/2016/california_ozone_plus_rcov_3.RDS")
## Add date covariates
epa <- addDoyField(epa)
epa <- addDowField(epa)
## Sites
sites <- getSites(epa)
## Scale target variable #####################################################
epa$sOzone <- scale(epa$Ozone)[,1]  ## scale returns a matrix
# epa$sOzone <- epa$sOzone + abs(min(epa$sOzone,na.rm=T))
# plot(epa$sOzone,type="l")


## Build design matrix ###############################################################
## Create additional columns for type factor (method B)
## The Rural factor is the base
## The other types are addicional information ???
epa$fac2 <- 0
epa$Suburban <- 0
epa$fac3 <- 0
epa$Urban <- 0
epa$fac2[epa$Location.Setting=="SUBURBAN"] <- 1
epa$Suburban[epa$Location.Setting=="SUBURBAN"] <- epa$Temperature[epa$Location.Setting=="SUBURBAN"]
epa$fac3[epa$Location.Setting=="URBAN"] <- 1
epa$Urban[epa$Location.Setting=="URBAN"] <- epa$Temperature[epa$Location.Setting=="URBAN"]
#aggregate(fac2~Location.Setting,epa,sum)
epa=epa[,names(epa) != "Location.Setting"]; # Del Location.Setting

covariates <- c("Temperature", "fac2", "Suburban", "fac3", "Urban")


## Split data #####################################################################
folds <- readRDS("output/folds.RDS")
## Fold(1)
epa.train <- epa[epa$Station.Code %in% sites$Station.Code[folds!=1],] 
epa.test <- epa[epa$Station.Code %in% sites$Station.Code[folds==1],]
epa.test$sOzoneH <- NA
#View(epa.train); View(epa.test)

## Fit model ########################################################################
n.samples <- 3000
n.burnin=1000
p = 6; # Number of parameters: MAIN+fac2+Suburban+fac3+Urban+(intercept)
beta.prior.mean <- rep(0, times=p)  # betas in 0 by default mean disabled betas???
beta.prior.precision <- matrix(100, nrow=p, ncol=p) # large precision means betas are not precise by default so they get narrow with sampling/fitting???
prior.shape <- -p/2  # 2; ?bayesLMConjugate
prior.rate <- 0; # 1; ?bayesLMConjugate
#fitting_file$obs = sqrt(fitting_file$obs); # SQRT transform 
## Notes: No transform!

## Notes:
## The obs variable of the fitting dataset was transformed 
## Obviously this is not necessary in the validation dataset
blm.1 <- bayesLMConjugate(sOzone~Temperature+fac2+Suburban+fac3+Urban, 
                          data=epa.train,
                          n.samples, beta.prior.mean, beta.prior.precision, 
                          prior.shape, prior.rate);
summary(blm.1$p.samples)
## 
## (WARNING) This is failing in the cluster
##
param = t(blm.1$p.samples[n.burnin:n.samples, 1:p]) #6x2001
sig.sq = blm.1$p.samples[n.burnin:n.samples, p+1]


## Calculate the fitted value = A*w ####################################
## In this case there are 2001 w's -> (param=6x2001)
## The first step:
## A*w generates (2001=samples-burnin) outputs per row (6 inputs)
## Then it is necesary to sample for each of those outputs
## using the appropiate distribution and SD!!!
if(forceRun){
  covars.fit = data.matrix(cbind(1,epa.train[,covariates]));
  fitted = covars.fit %*% param; # 9922 2001
  fitted_value = rep(0, nrow(fitted));
  ticToc({
    for(i in 1:ncol(fitted)){
      cat(".") #cat(i, "\n");
      sample = rnorm(n=nrow(fitted), mean=fitted[,i], sd=rep(sqrt(sig.sq[i]), nrow(fitted)));
      fitted_value = cbind(fitted_value, sample);
    }
    cat("\n")
  })
  fitted_value = fitted_value[,-1];
  saveRDS(fitted, file="output/spTimer/fitted.RDS")
  saveRDS(fitted_value, file="output/spTimer/fitted_value.RDS")
}else{
  fitted <- readRDS("output/spTimer/fitted.RDS")
  fitted_value <- readRDS("output/spTimer/fitted_value.RDS")
}

## Goodness of fit?
penalty = sum(apply(fitted_value,1,FUN=var)) # var per row
gft = sum((apply(fitted,1,FUN=mean)-epa.train$sOzone)^2,na.rm=T) #SSE
PMCC = gft + penalty;
## Notes
## Why one over fitted_value and the other over fitted?
## Maybe:
## fitted_value incorporates sig.sq -> estimates variability
## fitted is calculated using the mean of the parameters -> estimates mean output

## Validation? 
if(forceRun){
  covars.pred = data.matrix(cbind(1, epa.test[,covariates]))
  predicted = covars.pred %*% param
  predicted_value = rep(0, nrow(predicted))
  for(i in 1:ncol(predicted)){
    #cat(i, "\n");
    cat(".");
    sample = rnorm(n=nrow(predicted), mean=predicted[,i], sd=rep(sqrt(sig.sq[i]), nrow(predicted)));
    predicted_value = cbind(predicted_value, sample); # TRANSFORMATION!!!
  }
  predicted_value = predicted_value[,-1];
  saveRDS(predicted,file="output/spTimer/predicted.RDS")
  saveRDS(predicted_value,file="output/spTimer/predicted_value.RDS")
}else{
  predicted <- saveRDS("output/spTimer/predicted.RDS")
  predicted_value <- saveRDS("output/spTimer/predicted_value.RDS")
}

## Calculate confidence interval ############################################
lcl = apply(predicted_value, 1, quantile, 0.025); # lower confidence level
ucl = apply(predicted_value, 1, quantile, 0.975); # upper confidence level

# Count how many validation observation are inside the confidence interval
# This is called COVERAGE
total <- length(which(!is.na(epa.test$sOzone)));
count <- sum(lcl<=epa.test$sOzone & epa.test$sOzone<=ucl, na.rm=T)
prop = 100*count/total;  ## COVERAGE!

## Godness of fit 
## Average prediction
average_pred = apply(predicted_value, 1, FUN=mean)
evaluatePredictions(epa.test$sOzone, average_pred)
# diff =  epa.test$sOzone - average_pred
# rmse = sqrt(mean(diff^2, na.rm=T))
# mae = mean(abs(diff), na.rm=T)
# bias = mean(diff, na.rm=T);
# rbias = bias/mean(, na.rm=T);
plot(epa.test$sOzone, type="l")
lines(average_pred,col=2)

## Notes:
## Something is going bad!!!