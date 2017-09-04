## This code produces the parameter estimates of the best model 
## mentioned in Table 7. 
## The code requires installation of the software spAir. 

## Clean environment
rm(list=ls())

## Load libraries
library(spAir)


## Command line parameters ##########################################
useArgs <- T

if(useArgs){
  args = commandArgs(trailingOnly=TRUE)
  #print(args)
  
  if (length(args)<3) {
    stop("Three arguments must be supplied (nIter, nBurn, applyGLFilter?)", call.=FALSE)
  } else{
    nIter <- as.integer(args[1])
    nBurn <- as.integer(args[2])
    applyGLFilter<-as.logical(args[3])
  }
}


## Global variables #################################################
#applyGLFilter<-F
# iterations<- c(3000,5000,8000,10000)
# burns <- c(1000,1000,2000,2000)


###### reading the data files #######################################
spcheck1<-read.table("data/aurn/AURN_data_07_11.txt",header=T)


####
## Create a dataframe with the following columns:
## index, lon, lat, year, month, day, type, obs, sqrtaqm 
spcheck = spcheck1[,c("index","lon","lat","year","month","day","type")]
spcheck$obs=spcheck1$obs_no2 # the variable! 
spcheck$sqrtaqm = sqrt(spcheck1$aqum_no2) ## SQRT transform
## Order by index, date (year, month, day)
spcheck = spcheck[order(spcheck$index, spcheck$year, spcheck$month, spcheck$day),]


## Filter the data to reduce time period and spatial period ###################
if (applyGLFilter){
  spcheck <- spcheck[spcheck$lon>=-1.567 & spcheck$lon<=1.312
                     & spcheck$lat>=50.2 & spcheck$lat<=52.8,] #& spcheck$year==2011
  nrow(spcheck)
  length(unique(spcheck$index))
  ## Save DS V1
  # save(spcheck, file="data/aurn/AURN_data_GL_11.RData")
}


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
spcheck=spcheck[,-7]
#View(spcheck)


######## Choosing the fit and validation sites ########
siteIndexes <-unique(spcheck$index) 
sites <- unique(spcheck[,c("index","lon","lat")])

set.seed(11)

if(!applyGLFilter){
  # Original (tricky) method 
  sites.val <-c(4,9,28,30,32,44,50,75,102,108,111,122,130,175,194);
  sites.fit <- siteIndexes[!siteIndexes %in% sites.val]
  
  more.val.sites <- sample(sites.fit, size=(54-length(sites.val))) 
  sites.val <- c(sites.val, more.val.sites); 
  sites.fit <- siteIndexes[!siteIndexes %in% sites.val]
  
}else{
  # Alternative method
  sites.val <- sample(siteIndexes,15) #15
  sites.fit <- siteIndexes[!siteIndexes %in% sites.val] #47
}

## Quick scatter plot 
plot(lat~lon,sites[sites$index %in% sites.fit, ], col=1)
points(lat~lon,sites[sites$index %in% sites.val, ], col=2)


## Validation dataset #################################################################
prediction_file <- merge(data.frame(index=sites.val), spcheck)
prediction_file <- prediction_file[order(prediction_file$index,prediction_file$year,prediction_file$month,prediction_file$day),]
## Training dataset #################################################################
fitting_file<-merge(data.frame(index=sites.fit), spcheck)
fitting_file<-fitting_file[order(fitting_file$index,fitting_file$year,fitting_file$month,fitting_file$day),]
#fitting_file<-fitting_file[!is.na(fitting_file$obs),] # Drop NAs in obs
#View(prediction_file)
#View(fitting_file)


## All stations
coords_all = as.matrix(unique(spcheck[,2:3]))

## A grid ???
## How and why to generate this grid
## Idea, simply a grid at regular steps clipped over GB
## Prob??? A prior???
grid=read.table("data/aurn/grid_pop_eng.txt",header=F)
names(grid) <- c("lon","lat","prob")

## Quick plot
plot(lat~lon,grid, pch=4, col="gray")
points(lat~lon,sites[sites$index %in% sites.fit, ], col=2)
points(lat~lon,sites[sites$index %in% sites.val, ], col=3)

## Filter the grid too!
if(applyGLFilter){
  grid <- grid[grid$lon>=-1.567 & grid$lon<=1.312
                & grid$lat>=50.2 & grid$lat<=52.8,]
  ## Quick plot
  plot(lat~lon,grid, pch=4, col="gray")
  points(lat~lon,sites[sites$index %in% sites.fit, ], col=2)
  points(lat~lon,sites[sites$index %in% sites.val, ], col=3)
  ## Notes:
  ## The grid points are not regular!
}

## As matrix
cand_coords=as.matrix(grid[,1:2]) 


############################

cand_coords1=rep(0,(2*nrow(cand_coords)))
for(i in 1:nrow(cand_coords)){
  cand_coords1[(2*(i-1))+1] = cand_coords[i,1]  
  cand_coords1[(2*(i-1))+2] = cand_coords[i,2]
} 
grid_prob=grid[,3]
#grid_prob=rep(1, length(grid[,1])) # Alternate method???

## MCMC via Gibbs not using default choices
## hyper-parameters for the prior distributions
priors <- spT.priors(model="GPP",inv.var.prior=Gamm(2,1),beta.prior=Norm(0,10^10))
# initial values for the model parameters
initials <- spT.initials(model="GPP", sig2eps=0.01,sig2eta=0.5, beta=NULL, phi=0.001, sig2l=0.01)
# input for spatial decay
spatial.decay<-spT.decay(distribution="FIXED",value=0.001)

## Model fitting
set.seed(100)
coords.fit <- as.matrix(unique(fitting_file[,2:3]))
coords.val <- as.matrix(unique(prediction_file[,2:3]))
###########
# Define knots
coords.knots <- spT.grid.coords(Longitude=c(max(spcheck$lon),min(spcheck$lon)),
                                Latitude=c(max(spcheck$lat),min(spcheck$lat)), by=c(5,5))
plot(coords.knots, pch=3, col="gray")
points(lat~lon,sites[sites$index %in% sites.fit, ], col=2)
points(lat~lon,sites[sites$index %in% sites.val, ], col=3)

  posT <- c(0, 0)
  XY <- Formula.matrix(formula=obs ~ sqrtaqm + fac2 + Urban + fac3 + RKS, data=fitting_file)
  Y <- XY[[1]]
  X <- as.matrix(XY[[2]]) #TODO delete the as.matrix which is unnecessary
  time.data<-list(1,length(Y)/nrow(coords.fit))
  
  dates.fit = unique(fitting_file[,c(4:6)])
  #data_time$loc = 1:time.data[[2]]
  dates.fit$loc = 1:nrow(dates.fit)
  dates.val = unique(prediction_file[,c(4:6)])
  pos = which(dates.fit$year==dates.val$year[1] 
              & dates.fit$month==dates.val$month[1] 
              & dates.fit$day==dates.val$day[1])
  posT[1] = pos - 1

  end <- nrow(dates.val)
  pos <- which(dates.fit$year==dates.val$year[end] 
               & dates.fit$month==dates.val$month[end] 
               & dates.fit$day==dates.val$day[end])
  posT[2] = pos -1

### Fit model 
#for(i in 1:length(iterations)){
  # i <- 1
  #print(i)
  #nIter <- iterations[i]
  #nBurn <- burns[i]
  print("#########################################################")
  print("Fitting GPP model")
  model.output <- spT.Gibbs(formula=obs~sqrtaqm+fac2+Urban+fac3+RKS, data=fitting_file, model="GPP",
                            coords=coords.fit, knots.coords=coords.knots, newcoords=coords.val,
                            newdata=prediction_file, nItr=nIter, nBurn=nBurn,tol.dist=0.0001, #10000/5000
                            distance.method="geodetic:km",initials=initials,priors=priors,
                            scale.transform="SQRT",spatial.decay=spatial.decay,
                            cand_coords=cand_coords1,cand_len=nrow(cand_coords),
                            grid_prob=grid_prob,report=1000,predloc=posT)
  #print(model.output)
  areaName <- "UK"
  if(applyGLFilter){areaName <- "GL"}
  
  save(model.output,
       file=sprintf("output/model.output.%s.%d.RData",nIter))
  
  print("#########################################################")
  print("Results")
  spT.validation(prediction_file$obs,c(model.output$prediction[,1]))
#}

#if(!exists("model.output")){load("output/model.output.RData")}

# ## TODO
# load("output/model.output.UK.3000.RData")
# 
# 
# ##### produces the columns of table 7 ####################################
# summary(model.output)
# names(model.output)
# 
# parameter = model.output$parameters[,c(1,4,5)]
# rownames(parameter) = c("gamma_0", "gamma_1", "gamma_02", "gamma_12", "gamma_03", "gamma_13", "rho", "sigma_epsilon", "sigma_eta", "phi")
# parameter
# 
# 
# ## Metrics ###############################################################
# rmse=spT.validation(prediction_file$obs,c(model.output$prediction[,1]))


# #####  Producing the row of the Table 3 #############
# limits = read.table("data/aurn/Prediction_sites_summary.txt")
# limits = data.frame(mean=limits[,1], sd=limits[,2], low95=limits[,3], up95=limits[,4])
# limits = cbind(prediction_file[,c(1:7)], limits)
# limits = limits[which(!is.na(limits$obs)),]
# total_length = length(limits[,1])
# prop = 0
# for(ii in 1:total_length){
#   if((limits$low95[ii] < limits$obs[ii]) & (limits$obs[ii] < limits$up95[ii])){
#       prop = prop + 1
#   }
# }
# prop = 100*prop/total_length
# 
# ####### printing result will give row of Table 3 ############
# result = matrix(0, nrow=1, ncol=5)
# 
# result[1,1] = rmse[2];  result[1,2] = rmse[3]; result[1,3] = bias; result[1,4]=prop; result[1,5] = (cor(c(model.output$prediction[,1]), prediction_file$obs, use="pairwise.complete.obs"))^2;
# rownames(result) = c("Linear")
# colnames(result) = c("RMSPE", "MAPE", "Bias", "Coverage (%)", "R2")
# result



# Original script: 5000 / 1000