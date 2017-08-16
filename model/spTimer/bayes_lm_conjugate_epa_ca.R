## Paper [spTimer]
## ?bayesLMConjugate

## Clean the environment
rm(list=ls())

## Load libraries
library("spTimer")
library("akima")
library("coda")
library("spacetime")
library("fields")
library("forecast")
library("MASS")
library("mgcv")
library("spBayes")
library("colorspace") 
library("maps")
library("MBA")
library("openair")
source("util/my_helper.R")


## Read data #######################################################################
#epa <- readRDS("data/epa/epa_daily/2016/california_ozone.RDS")
epa <- readRDS("data/epa/epa_daily/2016/california_ozone_plus_rcov.RDS")
epa <- epa[order(epa$Station.Code, epa$Date),]
sites <- getSites(epa)
## Test
epa <- addJfield(epa)
epa$wd <- as.factor(weekdays(epa$Date))

## Standardize variable 
epa$sOzone <- scale(epa$Ozone)
epa$sOzone <- epa$sOzone + abs(min(epa$sOzone,na.rm=T))
plot(epa$sOzone,type="l")


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
folds <- readRDS("data/tmp/folds.RDS")
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

param = t(blm.1$p.samples[n.burnin:n.samples, 1:p]) #6x2001
sig.sq = blm.1$p.samples[n.burnin:n.samples, p+1]


## Calculate the fitted value = A*w ####################################
## In this case there are 2001 w's -> (param=6x2001)
## The first step:
## A*w generates (2001=samples-burnin) outputs per row (6 inputs)
## Then it is necesary to sample for each of those outputs
## using the appropiate distribution and SD!!!
if(T){
  covars.fit = data.matrix(cbind(1,epa.train[,covariates]));
  fitted = covars.fit %*% param; # 9922 2001
  fitted_value = rep(0, nrow(fitted));
  st <- Sys.time()
  for(i in 1:ncol(fitted)){
    cat(i, "\n");
    sample = rnorm(n=nrow(fitted), mean=fitted[,i], sd=rep(sqrt(sig.sq[i]), nrow(fitted)));
    fitted_value = cbind(fitted_value, sample);
  }
  Sys.time()-st
  fitted_value = fitted_value[,-1];
  save(fitted, fitted_value, file="model/spTimer/fitted.RData")
}
load("model/spTimer/fitted.RData")

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
if(T){
  covars.pred = data.matrix(cbind(1, epa.test[,covariates]))
  predicted = covars.pred %*% param
  predicted_value = rep(0, nrow(predicted))
  for(i in 1:ncol(predicted)){
    cat(i, "\n");
    sample = rnorm(n=nrow(predicted), mean=predicted[,i], sd=rep(sqrt(sig.sq[i]), nrow(predicted)));
    #predicted_value = cbind(predicted_value, sample^2); # TRANSFORMATION!!!
    predicted_value = cbind(predicted_value, sample); # TRANSFORMATION!!!
  }
  predicted_value = predicted_value[,-1];
  save(predicted,predicted_value,file="model/spTimer/predicted.RData")
}
load("model/spTimer/predicted.RData")

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