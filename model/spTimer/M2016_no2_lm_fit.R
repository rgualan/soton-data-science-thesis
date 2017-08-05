## Clean environment
rm(list=ls())

## Load libraries
library(spBayes)

## Reading the data files ######################################################
spcheck1<-read.table("data/aurn/AURN_data_07_11.txt",header=T)

####
## Create a dataframe with the following columns:
## index, lon, lat, year, month, day, type, obs, sqrtaqm 
spcheck = spcheck1[,c("index","lon","lat","year","month","day","type")];
spcheck$obs=spcheck1$obs_no2 # the variable! 
spcheck$sqrtaqm = sqrt(spcheck1$aqum_no2) ## SQRT transform
## Order by index, date (year, month, day)
spcheck = spcheck[order(spcheck$index, spcheck$year, spcheck$month, spcheck$day),];


## Filter the data to reduce time period and spatial period ###################
spcheck <- spcheck[spcheck$lon>=-1.567 & spcheck$lon<=1.312 
                    & spcheck$lat>=50.2 & spcheck$lat<=52.8
                    & spcheck$year==2011,]
nrow(spcheck)
length(unique(spcheck$index))
## Save DS V1
# save(spcheck, file="data/aurn/AURN_data_GL_11.RData")


## Preparing the data files ###################################################
## Create additional columns for type factor (method B)
spcheck$fac2 <- 0
spcheck$Urban <- 0
spcheck$fac3 <- 0
spcheck$RKS <- 0
spcheck$fac2[spcheck$type=="Urban"] <- 1
spcheck$Urban[spcheck$type=="Urban"] <- spcheck$sqrtaqm[spcheck$type=="Urban"]
spcheck$fac3[spcheck$type=="RKS"] <- 1
spcheck$RKS[spcheck$type=="RKS"] <- spcheck$sqrtaqm[spcheck$type=="RKS"]
spcheck=spcheck[,-7];
#View(spcheck)


######## Choosing the fit and validation sites ########
co_ind<-as.matrix(unique(spcheck[,1:3])); # index, lon, lat
nrow(co_ind) # 62

# Old method
# set.seed(11)
# #SitePred<-c(4,9,28,30,32,44,50,75,102,108,111,122,130,175,194); # predefined validation sites
# SitePred<-c(28,32,75,102,108,111,175); # predefined validation sites
# #sites for fitting
# gh<-co_ind[,1];
# sitefit<-gh[!gh %in% SitePred];
# new <- sample(sitefit, size=(54-length(SitePred)), replace = FALSE, prob = NULL); 
# SitePred <- append(SitePred, new); # 39 random additional validation sites
# sitefit<-gh[!gh %in% SitePred];
## Fit stations: 90
## Validation sites: 54 
## Proportion: 54/(90+54)

## New method
set.seed(11)
siteIndexes <-unique(spcheck$index) 
sites <- unique(spcheck[,c("index","lon","lat")])
sites.val <- sample(siteIndexes,15) #15
sites.fit <- siteIndexes[!siteIndexes %in% sites.val] #47
## Quick scatter plot 
plot(lat~lon,sites[sites$index %in% sites.fit, ], col=1)
points(lat~lon,sites[sites$index %in% sites.val, ], col=2)

## Validation dataset #################################################################
prediction_file <- merge(data.frame(index=sites.val), spcheck);
prediction_file <- prediction_file[order(prediction_file$index,prediction_file$year,prediction_file$month,prediction_file$day),];
## Training dataset #################################################################
fitting_file<-merge(data.frame(index=sites.fit), spcheck);
fitting_file<-fitting_file[order(fitting_file$index,fitting_file$year,fitting_file$month,fitting_file$day),];
fitting_file<-fitting_file[!is.na(fitting_file$obs),]; # Drop NAs in obs
#View(prediction_file)
#View(fitting_file)


## Fit model ########################################################################
n.samples <- 3000
n.burnin=1000
p = 6; # Number of parameters: sqrtaqm+fac2+Urban+fac3+RKS+(intercept)
beta.prior.mean <- rep(0, times=p)  # betas in 0 by default mean disabled betas???
beta.prior.precision <- matrix(100, nrow=p, ncol=p) # large precision means betas are not precise by default so they get narrow with sampling/fitting???
prior.shape <- 2; #??/
prior.rate <- 1; #??
fitting_file$obs = sqrt(fitting_file$obs); # SQRT transform


## Notes:
## The obs variable of the fitting dataset was transformed 
## Obviously this is not necessary in the validation dataset
blm.1 <- bayesLMConjugate(obs~sqrtaqm+fac2+Urban+fac3+RKS, data=fitting_file, 
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
covars.fit = data.matrix(cbind(1,fitting_file[,c(8:12)]));
fitted = covars.fit %*% param; # 9922 2001

# fitted_value = rep(0, nrow(fitted));
# st <- Sys.time()
# for(i in 1:ncol(fitted)){
#   cat(i, "\n");
#   sample = rnorm(n=nrow(fitted), mean=fitted[,i], sd=rep(sqrt(sig.sq[i]), nrow(fitted)));
#   fitted_value = cbind(fitted_value, sample);
# }
# Sys.time()-st
# fitted_value = fitted_value[,-1];
# save(fitted, fitted_value, file="data/tmp/fitted.RData")
load("data/tmp/fitted.RData")

## Goodness of fit?
penalty = sum(apply(fitted_value,1,FUN=var)) # var per row
gft = sum((apply(fitted,1,FUN=mean)-sqrt(fitting_file$obs))^2) #SSE
PMCC = gft + penalty;
## Notes
## Why one over fitted_value and the other over fitted?
## Maybe:
## fitted_value incorporates sig.sq -> estimates variability
## fitted is calculated using the mean of the parameters -> estimates mean output

## Validation? 
# covars.pred = data.matrix(cbind(1, prediction_file[,c(8:12)]))
# predicted = covars.pred %*% param
# predicted_value = rep(0, nrow(predicted))
# for(i in 1:ncol(predicted)){
#   cat(i, "\n");
#   sample = rnorm(n=nrow(predicted), mean=predicted[,i], sd=rep(sqrt(sig.sq[i]), nrow(predicted)));
#   predicted_value = cbind(predicted_value, sample^2); # TRANSFORMATION!!!
# }
# predicted_value = predicted_value[,-1];
# save(predicted,predicted_value,file="data/tmp/predicted.RData")
load("data/tmp/predicted.RData")


## Calculate confidence interval ############################################
lcl = apply(predicted_value, 1, quantile, 0.025); # lower confidence level
ucl = apply(predicted_value, 1, quantile, 0.975); # upper confidence level

# Count how many validation observation are inside the confidence interval
# This is called COVERAGE
total <- length(which(!is.na(prediction_file$obs)));
count <- sum(lcl<=prediction_file$obs & prediction_file$obs<=ucl, na.rm=T)
prop = 100*count/total;

## Godness of fit 
## Average prediction
average_pred = apply(predicted_value, 1, FUN=mean)
diff = average_pred - prediction_file$obs
rmse = sqrt(mean(diff^2, na.rm=T))
mae = mean(abs(diff), na.rm=T)
bias = mean(diff, na.rm=T);
rbias = bias/mean(prediction_file$obs, na.rm=T);

## Printing result will give row of Table 3 ############
result = matrix(0, nrow=1, ncol=5);
result[1,1] = rmse;  
result[1,2] = mae; 
result[1,3] = bias; 
result[1,4] = prop; 
result[1,5] = cor(average_pred, prediction_file$obs, use="pairwise.complete.obs")^2;
rownames(result) = c("Linear");
colnames(result) = c("RMSPE", "MAPE", "Bias", "Coverage (%)", "R2");
result


## Simple comparison versus a Linear Model
lm.1 = lm(obs~sqrtaqm+fac2+Urban+fac3+RKS, data=fitting_file) # collinearity???
#summary(lm.1)
diff.lm <- predict(lm.1, newdata = prediction_file)^2 - prediction_file$obs
rmse.lm = sqrt(mean(diff.lm^2, na.rm=T))
mae.lm = mean(abs(diff.lm), na.rm=T)
bias.lm = mean(diff.lm, na.rm=T)
rbias.lm = bias.lm/mean(prediction_file$obs, na.rm=T)
r2.lm = cor(predict(lm.1, newdata = prediction_file)^2, prediction_file$obs, use="pairwise.complete.obs")^2
result = rbind(result, lm=c(rmse.lm,mae.lm,bias.lm,NA,r2.lm))
result

## Notes:
## How is this possible???
## The Bayes LM performs a bit worst than lm
## despite its processing time

## Graphical comparison
plot(prediction_file$obs[1:100],type="l")
lines(predict(lm.1, newdata=prediction_file)^2, col=2)
lines(average_pred[1:100], col=3)
