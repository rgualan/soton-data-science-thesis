##### This code produces the parameter estimates of the best model mentioned in Table 7. The code requires installation of the software spAir. 

## Clean environment
rm(list=ls())

## Libraries ########################################################
library(spAir)
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
## Dates

## Sites
sites <- getSites(epa)
len_index=nrow(sites);
coords_all = as.matrix(sites[,c("UTM.X","UTM.Y")]);
index_data=sites$Station.Code;
len_index=length(index_data);

## Split data
folds <- getFolds()
## Fold(1)
epa.train <- epa[epa$Station.Code %in% sites$Station.Code[folds!=1],] 
epa.test <- epa[epa$Station.Code %in% sites$Station.Code[folds==1],]
epa.test$sOzoneH <- NA
#View(epa.train); View(epa.test)

## Specifying sample grid locations. One can just use his/her own choice. 
## Simple trick is to use   
knot.coords<-spT.grid.coords(
  Longitude=c(max(coords_all[,1]),min(coords_all[,1])),
  Latitude=c(max(coords_all[,2]),min(coords_all[,2])), 
  by=c(5,5))  
cand.coords<-spT.grid.coords(
  Longitude=c(max(coords_all[,1]),min(coords_all[,1])),
  Latitude=c(max(coords_all[,2]),min(coords_all[,2])), 
  by=c(25,25))  ## R should be a big number accroding to the instr
plot(knot.coords)
points(cand.coords, pch=".",col=2)
points(UTM.Y~UTM.X,sites,pch="+",col=3,cex=2)

cand_len = length(cand.coords[,1]);
cand_coords1=rep(0,(2*cand_len));
for(i in 1:cand_len){
  cand_coords1[(2*(i-1))+1] = cand.coords[i,1];  cand_coords1[(2*(i-1))+2] = cand.coords[i,2];
} 
grid_prob=rep(1, nrow(cand.coords));

## MCMC via Gibbs not using default choices
## hyper-parameters for the prior distributions
priors<-spT.priors(model="GPP",
                   inv.var.prior=Gamm(1,1),
                   beta.prior=Norm(0,10^10));
## initial values for the model parameters
initials<-spT.initials(model="GPP", 
                       sig2eps=1.01,
                       sig2eta=1.5, 
                       beta=NULL, 
                       phi=0.005, 
                       sig2l=0.01);
## input for spatial decay
spatial.decay<-spT.decay(distribution="FIXED",value=0.005);

####
set.seed(100)
coords.train<-unique(cbind(epa.train[,c("UTM.X","UTM.Y")]))
coords.test <-unique(cbind(epa.test[,c("UTM.X","UTM.Y")]))
###########
# Define knots
posT <- c(0, 0)
XY <- Formula.matrix(formula=sOzone ~ Temperature + fac2 + SuburbanTemp + fac3 + UrbanTemp, 
                     data=epa.train);
Y <- XY[[1]]
X <- as.matrix(XY[[2]])
time.data<-list(1,length(Y)/length(coords.train[,1]));
posT = c(0, time.data[[2]]-1)

## Train model
model.output <- spT.Gibbs(
  formula=sOzone ~ Temperature + fac2 + SuburbanTemp + fac3 + UrbanTemp,
  data=epa.train, model="GPP", coords=coords.train, knots.coords=knot.coords, 
  newcoords=coords.test, newdata=epa.test, nItr=10000, nBurn=5000,
  tol.dist=0.0001, distance.method="geodetic:km", initials=initials, 
  priors=priors, spatial.decay=spatial.decay, 
  cand_coords=cand_coords1, cand_len=cand_len, 
  grid_prob=grid_prob, report=1000, predloc=posT);
stop("debug!")

#scale.transform="SQRT",
saveRDS(model.output,file="data/model.output.RDS")
#readRDS("data/model.output.RDS")

rmse=spT.validation(prediction_file$obs_no2,c(model.output$prediction[,1]))
rmse
print(model.output)

##### produces the columns of table 7 ####################################
summary(model.output)
names(model.output)

parameter = model.output$parameters[,c(1,4,5)];
rownames(parameter) = c("gamma_0", "gamma_1", "gamma_02", "gamma_12", "gamma_03", "gamma_13", "rho", "sigma_epsilon", "sigma_eta", "phi");
parameter

#### produces the RMSPE and MAPE #########
rmse=spT.validation(prediction_file$obs_no2,c(model.output$prediction[,1]))
rmse

## Plot
# plot(prediction_file$obs_no2, type="l")
# lines(c(model.output$prediction[,1]), col=2)
