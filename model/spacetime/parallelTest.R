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
clusterExport(cl, c("timeLags","ticToc"))
#ticToc({
tryCatch({
  o <- clusterApply(cl, 1:length(timeLags), function(x){
    cat("Doing something in parallel: ", x, "\n")
    ticToc({
      for(i in 1:50) {
        a = rnorm(100000)
        a = sort(a)
      }
      print(length(a))
    })
    
    ## Step B
    try({
      print("Step B")
      stop("First problem")
    })
    
    ## Step C
    try({
      print("Step C")
      stop("Second problem")
    })

    out <- sum(1:x)
    #return(length(timeLags[[x]]))
    return(out)
  })
  print(o)
}, error = function(e){
  print(e)
}, finally = {
  stopCluster(cl)
})
#})

