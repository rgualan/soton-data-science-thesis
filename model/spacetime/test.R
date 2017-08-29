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
source("util/my_helper.R")

## Global variables
set.seed(123)
smplDays <- sort(sample(365,8))
forceRun <- T
paper <- setupPaper()
debugLevel <- T

## Read data ##############################################################################
epa <- readRDS("data/epa/epa_daily/2016/california_ozone_plus_rcov_3.RDS")
## Add date covariates
epa <- addDoyField(epa)
epa <- addDowField(epa)
## Sites
sites <- getSites(epa)
## Scale target variable #####################################################
epa$sOzone <- scale(epa$Ozone)[,1]  ## scale returns a matrix

## Assemble STFDF ############################################################
## Space dimension
epa.sp <- getSites(epa)
rownames(epa.sp) <- epa.sp$Station.Code
head(epa.sp)
coordinates(epa.sp) <- ~UTM.X+UTM.Y
proj4string(epa.sp) <- getUTMproj()

## Time dimension
epa.tm <- sort(unique(epa$Date))
epa.tm <- as.Date(epa.tm)  ## Ignore time data?? >> Corrects x labels problem in acf
# Combine the objects spatial, data-frame and time-dim into a STIDF:
epa.st <- STFDF(epa.sp,epa.tm,epa[,c("Date","Ozone","sOzone","Temperature",
                                     "RH","Rain","logRain","Wind","sqrtWind")]) 
# summary(epa.st)
# stplot(epa.st[,"2016-01-01::2016-01-08","sOzone"])
# dim(epa.st)


## Formula
#fm <- sOzone~1
fm <- sOzone ~ Temperature
fm <- sOzone ~ Temperature+RH+sqrtWind#+Elevation+Doy
#fm <- sOzone~Temperature+RH+Elevation

#a <- lm(fm,data=epa.st@data)


empVgm <- variogramST(fm, data=epa.st, na.omit=T)
(linStAni <- estiStAni(empVgm, c(10,500)))

prodSumModel <- vgmST("productSum",
                      space=vgm(0.03, "Exp", 75, 0.01),
                      time= vgm(0.2, "Sph", 7,  0.003), 
                      k=70)
printPlot(paper,"img/spacetime/vm_prodsum_eye.jpeg",7,5,FUN=function(){
  print(plot(empVgm, prodSumModel, all=T, map=F))
})
(fitProdSumModel <- fit.StVariogram(empVgm, prodSumModel, fit.method = 7, 
                                    stAni = linStAni, method = "L-BFGS-B", 
                                    control = list(parscale = c(1,10,1,1,0.1,1,10)),
                                    lower = rep(0.0001,7)))
attr(fitProdSumModel, "optim.output")$value
attr(fitProdSumModel, "MSE")
printPlot(paper,"img/spacetime/vm_prodsum_fit.jpeg",7,5,FUN=function(){
  print(plot(empVgm, fitProdSumModel, all=T, map=F))
})
# Exp+Exp: 0.002605351
# Exp+Sph: 
# Sph+Exp: 
# Sph+Sph: 
printPlot(paper,"img/spacetime/vm_prodsum_fit_2.jpeg",5,5,FUN=function(){
  print(plot(empVgm, fitProdSumModel, wireframe=T, all=T, scales=list(arrows=F),zlab="'"))
})




## 10-fold cross-validation ######################################################################
runFold <- function(k, folds, fm, epa.st, outputColumnName, model, linStAni, nmax=50){
  res <- matrix(NA, dim(epa.st)[1], dim(epa.st)[2])
  ticToc({
    for(i in 1:k) {
      cat("Fold", i, "\n")
      testIndices <- folds==i
      
      for(j in which(folds==i)){
        tmp <- krigeST(fm, data=epa.st[!testIndices],
                       newdata=epa.st[j,drop=F], 
                       model, nmax=nmax,
                       stAni=linStAni/24/3600,
                       progress=F)$var1.pred
        res[j, !is.na(epa.st[j,])[,"sOzone"] ] <- tmp
        # Simple partial assessment
        if (debugLevel) print(evaluatePredictions(as.vector(epa.st[j,,"sOzone"][,1]), res[j,]))
        printPlot(paper,paste0("img/spacetime/cv/site",j,".jpg"),6,5,FUN=function(){
          plot(as.numeric(epa.st[j,,"sOzone"][,1]), type="l",ylab=""); lines(res[j,],col=2)
        })
        # cat("=")
      }
      # cat("\n")
      #stop("Debug!")
    }
    epa.st@data$yh <- as.vector(res)
    names(epa.st@data)[ncol(epa.st@data)] <- outputColumnName
  })
  return(epa.st)
}

k <- 1 # TODO!
folds <- readRDS("output/folds.RDS")
epa.st <- runFold(k, folds, fm, epa.st, "psModel", fitProdSumModel, linStAni, nmax=50)
#saveRDS(epa.st,file="output/spacetime/epa.st.RDS")

## CV evaluation metrics
evaluatePredictions(epa.st[,,"psModel",drop=F]@data[[1]], 
                      epa.st[,,"sOzone",drop=F]@data[[1]])

