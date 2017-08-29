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
forceRun <- T # Run heavy weight processing
paper <- setupPaper() # Print the plots in image format for the paper
debugLevel <- T # Present additional debugging information

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
epa.st <- assembleSTFDF(epa) 
# stplot(epa.st[,"2016-01-01::2016-01-08","sOzone"])
# dim(epa.st)

## Run simulations in parallel ###############################################################
## Objective:
## Try several covariates 
## To assess which combination produces the best results
## For time restrictions, the assessment is made only on the first fold
## Several formulas
## Each formula trying three sets of time lags
## Each experiment tries the four covariance models
formulas <- c("sOzone~Elevation",
              "sOzone~Elevation+Temperature",
              "sOzone~Elevation+Temperature+RH",
              "sOzone~Elevation+Temperature+RH+sqrtWind")
timeLags <- list(a=0:5, b=0:7, c=0:10)

cl <- makeCluster(4, outfile="")
clusterExport(cl, c("epa.st","paper","timeLags","formulas"))
tryCatch({
  clusterApply(cl, 1:length(formulas), function(idx){
    source("model/spacetime/testStKriging.R")
    print(idx)
    try({
      m <- test_st_kriging_using_1fold(formulas[idx], epa.st, timeLags[[1]], paper, as.character(idx))
      print(sprintf("Fm: %s. Tl: %d",formulas[idx],1))
      print(m)
    })
    try({
      m <- test_st_kriging_using_1fold(formulas[idx], epa.st, timeLags[[2]], paper, as.character(idx))
      print(sprintf("Fm: %s. Tl: %d",formulas[idx],2))
      print(m)
    })
    try({
      m <- test_st_kriging_using_1fold(formulas[idx], epa.st, timeLags[[3]], paper, as.character(idx))
      print(sprintf("Fm: %s. Tl: %d",formulas[idx],3))
      print(m)
    })
  })
}, error = function(e){
  print(e)
}, finally = {
  stopCluster(cl)
})
