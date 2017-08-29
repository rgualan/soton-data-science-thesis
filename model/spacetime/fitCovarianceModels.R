
## Function that fits the 4 simplest Covariance models
## implemented in gstat
fitCovarianceModels <- function(epa.st, fm, tlags, paper){
  ## WORKFLOW ##################################################################
  empVgm <- variogramST(fm, data=epa.st, tlags=tlags, na.omit=T)
  linStAni <- estiStAni(empVgm, c(10,200))
  printPlot(paper,"img/spacetime/empirical_variogram.jpeg",6,6,FUN=function(){
    print(plot(empVgm, map=F))
  })

  
  ## Separable #################################################################
  separableModel <- vgmST("separable",
                          space = vgm(0.45,"Exp",500,0.5), ##Factor,function,range,nugget
                          time = vgm(0.45,"Exp",7,0.5), ##Factor,function,range,nugget 
                          sill=0.8) # Base
  printPlot(paper,"img/spacetime/vm_separable_eye.jpeg",7,5,FUN=function(){
    print(plot(empVgm, separableModel, all=T, map=F, zlab=""))
  })
  (fitSepModel <- fit.StVariogram(empVgm, separableModel))
  attr(fitSepModel, "optim.output")$value
  attr(fitSepModel, "MSE")
  printPlot(paper,"img/spacetime/vm_separable_fit.jpeg",7,5,FUN=function(){
    print(plot(empVgm, fitSepModel, all=T, map=F))
  })
  printPlot(paper,"img/spacetime/vm_separable_fit_2.jpeg",5,5,FUN=function(){
    print(plot(empVgm, fitSepModel, wireframe=T, all=T, scales=list(arrows=F), zlab=""))
  })
  
  ## Product-sum ###############################################################
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
  printPlot(paper,"img/spacetime/vm_prodsum_fit_2.jpeg",5,5,FUN=function(){
    print(plot(empVgm, fitProdSumModel, wireframe=T, all=T, scales=list(arrows=F),zlab="'"))
  })
  
  ## Metric ################################################################################
  metricModel <- vgmST("metric",
                       joint=vgm(0.4, "Mat", 100, 0.3, kappa = 1.5),
                       stAni=50)
  printPlot(paper,"img/spacetime/vm_metric_eye.jpeg",7,5,FUN=function(){
    print(plot(empVgm, metricModel, all=T, map=F))
  })
  (fitMetricModel <- fit.StVariogram(empVgm, metricModel, fit.method = 7,
                                     stAni = linStAni, method = "L-BFGS-B"))
  attr(fitProdSumModel, "optim.output")$value
  attr(fitProdSumModel, "MSE")
  printPlot(paper,"img/spacetime/vm_metric_fit.jpeg",7,5,FUN=function(){
    print(plot(empVgm, fitProdSumModel, all=T, map=F))
  })
  printPlot(paper,"img/spacetime/vm_metric_fit_2.jpeg",5,5,FUN=function(){
    print(plot(empVgm, fitProdSumModel, wireframe=T, all=T, scales=list(arrows=F),zlab=""))
  })
  
  ## Sum-metric ############################################################################
  sumMetricModel <- vgmST("sumMetric",
                          space = vgm(0.15, "Sph", 75, 0.07),
                          time = vgm(0.16, "Exp", 7, 0.1),
                          joint = vgm(0.25, "Sph", 150, 0.1),
                          stAni = linStAni)
  printPlot(paper,"img/spacetime/vm_summetric_eye.jpeg",7,5,FUN=function(){
    print(plot(empVgm, sumMetricModel, all=T, map=F))
  })
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
  printPlot(paper,"img/spacetime/vm_summetric_fit.jpeg",7,5,FUN=function(){
    print(plot(empVgm, fitSumMetricModel, all=T, map=F))
  })
  printPlot(paper,"img/spacetime/vm_summetric_fit_2.jpeg",5,5,FUN=function(){
    print(plot(empVgm, fitSumMetricModel, wireframe=T, all=T, scales=list(arrows=F),zlab=""))
  })
  
  ## All the models in one Figure ###################################################
  # A wireframe (3D) plot of sample variogram and fitted variogram models can
  # be obtained e.g. by
  printPlot(paper,"img/spacetime/all.jpeg",7,5,FUN=function(){
    print(plot(vv, list(fitSepModel, fitProdSumModel, fitMetricModel, fitSumMetricModel), 
               all=T, wireframe=T, #zlim=c(0,120),
               zlab=NULL,
               xlab=list("Distance (km)", rot=30, cex=0.6),
               ylab=list("Time lag (days)", rot=-35, cex=0.6),
               scales=list(arrows=F, z = list(distance = 5), cex=0.5)))
  })
  

  ## Assemble output
  out <- list(separable=fitSepModel,
              prodSum=fitProdSumModel,
              metric=fitMetricModel,
              sumMetric=fitSumMetricModel,
              linStAni=linStAni)
  
  return(out)  
}
