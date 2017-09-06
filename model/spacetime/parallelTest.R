## Libraries
library(gstat)
library(spacetime)
library(reshape2)
library(parallel)
source("util/my_helper.R")
source("model/spacetime/testStKriging.R")


## Run simulations in parallel
cl <- makeCluster(3,outfile="")
timeLags <- list(a=0:5, b=0:7, c=0:10)
mainDf <- data.frame(id=1:10)
clusterExport(cl, c("timeLags","ticToc","mainDf"))
tryCatch({
  o <- clusterApply(cl, 1:length(timeLags), function(x){
    cat("Doing something in parallel: ", x, "\n")
    mainDf[sprintf("column_%d",x)] <- mainDf$id*x
    return(mainDf)
  })
}, error = function(e){
  print(e)
}, finally = {
  stopCluster(cl)
})

out = Reduce(function(...) merge(..., all=T), o)
head(out)
