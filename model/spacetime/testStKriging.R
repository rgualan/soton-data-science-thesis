## Ref: 
## Spatio-Temporal data in R
## https://www.r-bloggers.com/spatio-temporal-kriging-in-r/
## https://cran.r-project.org/web/packages/gstat/vignettes/st.pdf

## Libraries
source("util/my_helper.R")
library(gstat)
library(spacetime)
library(reshape2)

## Global variables
debugLevel=F


## 10-fold cross-validation ######################################################################
runCV <- function(k, folds, fm, epa.st, outputColumnName, model, linStAni, nmax=50, expName=NA){
  res <- matrix(NA, dim(epa.st)[1], dim(epa.st)[2])
  ticToc({
    for(i in 1:k) { #1:k
      cat("Fold", i, "\n")
      testIndices <- folds==i
      
      for(j in which(folds==i)){
        #cat("Station.Code: ",as.character(epa.st@sp[j,]$Station.Code), "\n")
        tmp <- krigeST(fm, data=epa.st[!testIndices],
                       newdata=epa.st[j,drop=F], 
                       model, nmax=nmax,
                       stAni=linStAni/24/3600,
                       progress=F)$var1.pred
        res[j, !is.na(epa.st[j,])[,"sOzone"] ] <- tmp
        # Simple partial assessment
        if (debugLevel){
          print(evaluatePredictions(as.vector(epa.st[j,,"sOzone"][,1]), res[j,]))
          tmpFolder <- "img/spacetime/cv/"
          if(!is.na(expName)){
            #if( !dir.exists(paste0(tmpFolder,expName)) ){
            if( !file.exists(paste0(tmpFolder,expName)) ){
              dir.create(file.path(tmpFolder,expName))
            }
            tmpFolder <- paste0(tmpFolder, expName, "/")
          }
          printPlot(paper,paste0(tmpFolder,as.character(epa.st@sp[j,]$Station.Code),".jpg"),6,5,FUN=function(){
            plot(as.numeric(epa.st[j,,"sOzone"][,1]), type="l",ylab=""); lines(res[j,],col=2)
          })
        } 
        # cat("=")
      }
      # cat("\n")
      # stop("Debug!")
    }
    epa.st@data$yh <- as.vector(res)
    names(epa.st@data)[ncol(epa.st@data)] <- outputColumnName
  })
  return(epa.st)
}


## Complete st workflow #####################################################################
test_st_kriging_using_10cv <- function(fm, epa.st, tlags, paper, expName=NA){
  empVgm <- variogramST(fm, data=epa.st, tlags=tlags, na.omit=T)

  linStAni <- estiStAni(empVgm, c(10,200))
  
  ## Separable
  separableModel <- vgmST("separable",
                          space = vgm(0.45,"Exp",500,0.5), ##Factor,function,range,nugget
                          time = vgm(0.45,"Exp",7,0.5), ##Factor,function,range,nugget 
                          sill=0.8) # Base
  (fitSepModel <- fit.StVariogram(empVgm, separableModel))
  attr(fitSepModel, "optim.output")$value
  attr(fitSepModel, "MSE")

  ## Product-sum
  prodSumModel <- vgmST("productSum",
                        space=vgm(0.03, "Exp", 75, 0.01),
                        time= vgm(0.2, "Sph", 7,  0.003), 
                        k=70)
  (fitProdSumModel <- fit.StVariogram(empVgm, prodSumModel, fit.method = 7, 
                                      stAni = linStAni, method = "L-BFGS-B", 
                                      control = list(parscale = c(1,10,1,1,0.1,1,10)),
                                      lower = rep(0.0001,7)))
  attr(fitProdSumModel, "optim.output")$value
  attr(fitProdSumModel, "MSE")

  ## Sum-metric
  sumMetricModel <- vgmST("sumMetric",
                          space = vgm(0.15, "Sph", 75, 0.07),
                          time = vgm(0.16, "Exp", 7, 0.1),
                          joint = vgm(0.25, "Sph", 150, 0.1),
                          stAni = linStAni)
  (fitSumMetricModel <- 
      fit.StVariogram(empVgm, sumMetricModel, fit.method = 7, stAni=linStAni,
                      method = "L-BFGS-B", 
                      lower = c(sill.s = 0,  range.s = 10,  nugget.s = 0,
                                sill.t = 0,  range.t = 0.1, nugget.t = 0,
                                sill.st= 0,  range.st = 10, nugget.st = 0, 
                                anis = 40),
                      upper = c(sill.s = 10, range.s = 500,  nugget.s = 5,
                                sill.t = 10, range.t = 15,   nugget.t = 5,
                                sill.st= 10, range.st = 500, nugget.st = 5,
                                anis = 500),
                      control = list(parscale = c(1,100,1,1,0.5,1,1,100,1,100),
                                     maxit=1e5)))
  attr(fitSumMetricModel, "optim.output")$value
  attr(fitSumMetricModel, "MSE")

  ## 10-fold CV  
  k <- 10
  folds <- getFolds()
  epa.st <- runFold(k, folds, fm, epa.st, "sepModel", fitSepModel, linStAni, nmax=50, expName)
  epa.st <- runFold(k, folds, fm, epa.st, "psModel", fitProdSumModel, linStAni, nmax=50, expName)
  epa.st <- runFold(k, folds, fm, epa.st, "sumMetricModel", fitSumMetricModel, linStAni, nmax=50, expName) ## TODO: the var
  epa.st <- runFold(k, folds, fm, epa.st, "prodSumModel", fitProdSumModel, linStAni, nmax=50, expName)

  ## CV evaluation metrics
  rbind(
    evaluatePredictions(epa.st[,,"sepModel",drop=F]@data[[1]], 
                        epa.st[,,"sOzone",drop=F]@data[[1]]),
    evaluatePredictions(epa.st[,,"psModel",drop=F]@data[[1]], 
                        epa.st[,,"sOzone",drop=F]@data[[1]]),
    evaluatePredictions(epa.st[,,"sumMetricModel",drop=F]@data[[1]], 
                        epa.st[,,"sOzone",drop=F]@data[[1]]),
    evaluatePredictions(epa.st[,,"prodSumModel",drop=F]@data[[1]], 
                        epa.st[,,"sOzone",drop=F]@data[[1]])
  )

  return(epa.st)
}


## Single fold - workflow #####################################################################
test_st_kriging_using_1fold <- function(fm, epa.st, tlags, paper, expName=NA){
  if(class(fm)!="formula"){
    fm <- as.formula(fm) 
  }
  empVgm <- variogramST(fm, data=epa.st, tlags=tlags, na.omit=T)
  linStAni <- estiStAni(empVgm, c(10,200))
  
  ## Separable
  separableModel <- vgmST("separable",
                          space = vgm(0.45,"Exp",500,0.5), ##Factor,function,range,nugget
                          time = vgm(0.45,"Exp",7,0.5), ##Factor,function,range,nugget 
                          sill=0.8) # Base
  (fitSepModel <- fit.StVariogram(empVgm, separableModel))
  attr(fitSepModel, "optim.output")$value
  attr(fitSepModel, "MSE")
  
  ## Product-sum
  prodSumModel <- vgmST("productSum",
                        space=vgm(0.03, "Exp", 75, 0.01),
                        time= vgm(0.2, "Sph", 7,  0.003), 
                        k=70)
  (fitProdSumModel <- fit.StVariogram(empVgm, prodSumModel, fit.method = 7, 
                                      stAni = linStAni, method = "L-BFGS-B", 
                                      control = list(parscale = c(1,10,1,1,0.1,1,10)),
                                      lower = rep(0.0001,7)))
  attr(fitProdSumModel, "optim.output")$value
  attr(fitProdSumModel, "MSE")
  
  ## Metric
  metricModel <- vgmST("metric",
                       joint=vgm(0.4, "Mat", 100, 0.3, kappa = 1.5),
                       stAni=50)
  (fitMetricModel <- fit.StVariogram(empVgm, metricModel, fit.method = 7,
                                     stAni = linStAni, method = "L-BFGS-B"))

  ## Sum-metric
  sumMetricModel <- vgmST("sumMetric",
                          space = vgm(0.15, "Sph", 75, 0.07),
                          time = vgm(0.16, "Exp", 7, 0.1),
                          joint = vgm(0.25, "Sph", 150, 0.1),
                          stAni = linStAni)
  (fitSumMetricModel <- 
      fit.StVariogram(empVgm, sumMetricModel, fit.method = 7, stAni=linStAni,
                      method = "L-BFGS-B", 
                      lower = c(sill.s = 0,  range.s = 10,  nugget.s = 0,
                                sill.t = 0,  range.t = 0.1, nugget.t = 0,
                                sill.st= 0,  range.st = 10, nugget.st = 0, 
                                anis = 40),
                      upper = c(sill.s = 10, range.s = 500,  nugget.s = 5,
                                sill.t = 10, range.t = 15,   nugget.t = 5,
                                sill.st= 10, range.st = 500, nugget.st = 5,
                                anis = 500),
                      control = list(parscale = c(1,100,1,1,0.5,1,1,100,1,100),
                                     maxit=1e5)))
  attr(fitSumMetricModel, "optim.output")$value
  attr(fitSumMetricModel, "MSE")
  
  ## 10-fold CV  
  k <- 10
  folds <- getFolds()
  epa.st <- runFold(1, folds, fm, epa.st, "sepModel", fitSepModel, linStAni, nmax=50, expName)
  epa.st <- runFold(1, folds, fm, epa.st, "metricModel", fitMetricModel, linStAni, nmax=50, expName)
  epa.st <- runFold(1, folds, fm, epa.st, "psModel", fitProdSumModel, linStAni, nmax=50, expName)
  epa.st <- runFold(1, folds, fm, epa.st, "sumMetricModel", fitSumMetricModel, linStAni, nmax=50, expName) ## TODO: the var
  
  
  ## CV evaluation metrics
  metrics <- rbind(
    evaluatePredictions(epa.st[,,"sepModel",drop=F]@data[[1]], 
                        epa.st[,,"sOzone",drop=F]@data[[1]]),
    evaluatePredictions(epa.st[,,"metricModel",drop=F]@data[[1]],
                        epa.st[,,"sOzone",drop=F]@data[[1]]),
    evaluatePredictions(epa.st[,,"psModel",drop=F]@data[[1]], 
                        epa.st[,,"sOzone",drop=F]@data[[1]]),
    evaluatePredictions(epa.st[,,"sumMetricModel",drop=F]@data[[1]], 
                        epa.st[,,"sOzone",drop=F]@data[[1]])
  )
  
  return(metrics)
}

