library(xtable)

## Function that fits the 4 simplest Covariance models
## implemented in gstat
fitCovarianceModels <- function(epa.st, fm, tlags, paper){
  ## WORKFLOW ##################################################################
  empVgm <- variogramST(fm, data=epa.st, tlags=c(0:2), na.omit=T)
  linStAni <- estiStAni(empVgm, c(10,200))
  printPlot(paper,"img/spacetime/empirical_variogram.jpeg",6,6,FUN=function(){
    print(plot(empVgm, map=F))
  })

  
  ## Separable #################################################################
  print("Separable")
  separableModel <- vgmST("separable",
                          space = vgm(0.7,"Exp",100,0.1), ##Factor,function,range,nugget
                          time = vgm(0.7,"Exp",10,0.5), ##Factor,function,range,nugget 
                          sill=0.3) # Base

  printPlot(paper,"img/spacetime/vm_separable_eye.jpeg",7,5,FUN=function(){
    print(plot(empVgm, separableModel, all=T, map=F, zlab=""))
  })
  (fitSepModel <- fit.StVariogram(empVgm, separableModel))
  attr(fitSepModel, "optim.output")$value
  print(format(attr(fitSepModel, "MSE"), scientific=T, digits = 3))
  printPlot(paper,"img/spacetime/vm_separable_fit.jpeg",7,5,FUN=function(){
    print(plot(empVgm, fitSepModel, all=T, map=F))
  })
  printPlot(paper,"img/spacetime/vm_separable_fit_2.jpeg",5,5,FUN=function(){
    print(plot(empVgm, fitSepModel, wireframe=T, all=T, scales=list(arrows=F), zlab=""))
  })
  ## Summary table for the report
  a<-fitSepModel[[1]]; b<-fitSepModel[[2]]
  sumTable <- 
    rbind(space=data.frame("psill"=a[1,2], model=a[2,1], range=a[2,3], nugget=a[1,3], stsill=fitSepModel[[3]]),
          time=data.frame("psill"=b[1,2], model=b[2,1], range=b[2,3], nugget=b[1,3], stsill=NA))
  print(sumTable); print(xtable(sumTable))
  
  ## Product-sum ###############################################################
  print("Product-sum")
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
  print(format(attr(fitProdSumModel, "MSE"), scientific=T, digits = 3))
  printPlot(paper,"img/spacetime/vm_prodsum_fit.jpeg",7,5,FUN=function(){
    print(plot(empVgm, fitProdSumModel, all=T, map=F))
  })
  printPlot(paper,"img/spacetime/vm_prodsum_fit_2.jpeg",5,5,FUN=function(){
    print(plot(empVgm, fitProdSumModel, wireframe=T, all=T, scales=list(arrows=F),zlab="'"))
  })
  a<-fitProdSumModel[[1]]; b<-fitProdSumModel[[2]]
  sumTable <- 
    rbind(space=data.frame("psill"=a[1,2], model=a[2,1], range=a[2,3], nugget=a[1,3], k=fitSepModel[[3]]),
          time=data.frame("psill"=b[1,2], model=b[2,1], range=b[2,3], nugget=b[1,3], k=NA))
  print(sumTable); print(xtable(sumTable))
  
  ## Metric ################################################################################
  print("Metric")
  metricModel <- vgmST("metric",
                       joint=vgm(0.2, "Exp", 300, 0.1),
                       stAni=200)
  printPlot(paper,"img/spacetime/vm_metric_eye.jpeg",7,5,FUN=function(){
    print(plot(empVgm, metricModel, all=T, map=F))
  })
  (fitMetricModel <- fit.StVariogram(empVgm, metricModel, fit.method = 7,
                                     stAni = linStAni, method = "L-BFGS-B"))
  attr(fitMetricModel, "optim.output")$value
  print(format(attr(fitMetricModel, "MSE"), scientific=T, digits = 3))
  printPlot(paper,"img/spacetime/vm_metric_fit.jpeg",7,5,FUN=function(){
    print(plot(empVgm, fitMetricModel, all=T, map=F))
  })
  printPlot(paper,"img/spacetime/vm_metric_fit_2.jpeg",5,5,FUN=function(){
    print(plot(empVgm, fitMetricModel, wireframe=T, all=T, scales=list(arrows=F),zlab=""))
  })
  ## Summary table for the report
  a<-fitMetricModel[[1]];
  sumTable <- 
    rbind(joint=data.frame("psill"=a[1,2], model=a[2,1], range=a[2,3], nugget=a[1,3], anisotropy=fitMetricModel[[2]]))
  print(sumTable); print(xtable(sumTable))
  
  ## Sum-metric ############################################################################
  print("Sum-metric")
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
  print(format(attr(fitSumMetricModel, "MSE"), scientific=T, digits = 3))
  printPlot(paper,"img/spacetime/vm_summetric_fit.jpeg",7,5,FUN=function(){
    print(plot(empVgm, fitSumMetricModel, all=T, map=F))
  })
  printPlot(paper,"img/spacetime/vm_summetric_fit_2.jpeg",5,5,FUN=function(){
    print(plot(empVgm, fitSumMetricModel, wireframe=T, all=T, scales=list(arrows=F),zlab=""))
  })
  ## Summary table for the report
  a<-fitSumMetricModel[[1]]; b<-fitSumMetricModel[[2]]; c<-fitSumMetricModel[[3]];
  sumTable <- 
    rbind(space=data.frame("psill"=a[1,2], model=a[2,1], range=a[2,3], nugget=a[1,3], anisotropy=fitSumMetricModel[[4]]),
          time=data.frame("psill"=b[1,2], model=b[2,1], range=b[2,3], nugget=b[1,3], anisotropy=NA),
          joint=data.frame("psill"=c[1,2], model=c[2,1], range=c[2,3], nugget=c[1,3], anisotropy=NA))
  print(sumTable); print(xtable(sumTable))
  
  ## All the models in one Figure ###################################################
  # A wireframe (3D) plot of sample variogram and fitted variogram models can
  # be obtained e.g. by
  printPlot(paper,"img/spacetime/all.jpeg",7,5,FUN=function(){
    print(plot(empVgm, list(fitSepModel, fitProdSumModel, fitMetricModel, fitSumMetricModel), 
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
