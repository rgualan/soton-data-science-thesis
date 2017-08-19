## REF: http://s3l.stanford.edu/blog/?p=73

library(gbm)

set.seed(100)
N <- 1000
X1 <- runif(N)
X2 <- 2*runif(N)
X3 <- ordered(sample(letters[1:4],N,replace=TRUE),levels=letters[4:1])
X4 <- factor(sample(letters[1:6],N,replace=TRUE))
X5 <- factor(sample(letters[1:3],N,replace=TRUE))
X6 <- 3*runif(N)
mu <- c(-1,0,1,2)[as.numeric(X3)]

SNR <- 10 # signal-to-noise ratio
Y <- X1**1.5 + 2 * (X2**.5) + mu
sigma <- sqrt(var(Y)/SNR)
Y <- Y + rnorm(N,0,sigma)

data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)
head(data)



gbm1 <- gbm(Y~X1+X2+X3+X4+X5+X6,         # formula
            data=data,                   # dataset
            var.monotone=c(0,0,0,0,0,0), # -1: monotone decrease,
            # +1: monotone increase,
            #  0: no monotone restrictions
            distribution="gaussian",     # bernoulli, adaboost, gaussian,
            # poisson, coxph, and quantile available
            n.trees=3000,                # number of trees
            shrinkage=0.005,             # shrinkage or learning rate,
            # 0.001 to 0.1 usually work
            interaction.depth=3,         # 1: additive model, 2: two-way interactions, etc.
            bag.fraction = 0.5,          # subsampling fraction, 0.5 is probably best
            train.fraction = 0.5,        # fraction of data for training,
            # first train.fraction*N used for training
            n.minobsinnode = 10,         # minimum total weight needed in each node
            cv.folds = 5,                # do 5-fold cross-validation
            keep.data=TRUE,              # keep a copy of the dataset with the object
            verbose=TRUE)                # print out progress

# check performance using an out-of-bag estimator
# OOB underestimates the optimal number of iterations
best.iter <- gbm.perf(gbm1,method="OOB")
print(best.iter)

data.predict = predict(gbm1, n.trees = best.iter)
SSE = sum((Y - data.predict)^2)
print(SSE)


# =========================================================================
X1[sample(1:N,size=500)] <- NA
X4[sample(1:N,size=300)] <- NA

data2 <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)

gbm2 <- gbm(Y~X1+X2+X3+X4+X5+X6,         # formula
            data=data2,                   # dataset
            var.monotone=c(0,0,0,0,0,0), # -1: monotone decrease,
            # +1: monotone increase,
            #  0: no monotone restrictions
            distribution="gaussian",     # bernoulli, adaboost, gaussian,
            # poisson, coxph, and quantile available
            n.trees=3000,                # number of trees
            shrinkage=0.005,             # shrinkage or learning rate,
            # 0.001 to 0.1 usually work
            interaction.depth=3,         # 1: additive model, 2: two-way interactions, etc.
            bag.fraction = 0.5,          # subsampling fraction, 0.5 is probably best
            train.fraction = 0.5,        # fraction of data for training,
            # first train.fraction*N used for training
            n.minobsinnode = 10,         # minimum total weight needed in each node
            cv.folds = 5,                # do 5-fold cross-validation
            keep.data=TRUE,              # keep a copy of the dataset with the object
            verbose=TRUE)                # print out progress

# check performance using an out-of-bag estimator
# OOB underestimates the optimal number of iterations
best.iter <- gbm.perf(gbm2,method="OOB")
print(best.iter)

data.predict = predict(gbm2, n.trees = best.iter)
SSE = sum((Y - data.predict)^2)
print(SSE)


# =============================================================================
gbmImpute = function(x, max.iters = 2, cv.fold = 5, verbose=T) {
  
  missing.matrix = is.na(x)
  numMissing = sum(missing.matrix)
  if(verbose) {
    print(paste("imputing on", numMissing, "missing values with matrix size",
                nrow(x)*ncol(x), sep=" "))
  }
  if(numMissing == 0) {
    return (x)
  }
  
  missing.cols.indices = which(apply(missing.matrix, 2, function(i) {
    any(i)
  }))
  
  for (i in 1:max.iters) {
    if (verbose) print (paste("Begin iteration: ", i))
    x[,missing.cols.indices] = sapply(missing.cols.indices, function(j) {
      good.data = which(!missing.matrix[,j])
      gbm1 <- gbm(x[good.data,j] ~ .,
                  data = as.data.frame(x[good.data,-j]),
                  var.monotone = rep(0, ncol(x)-1), # -1: monotone decrease,
                  # +1: monotone increase,
                  #  0: no monotone restrictions
                  distribution="gaussian",     # bernoulli, adaboost, gaussian,
                  # poisson, coxph, and quantile available
                  n.trees=3000,                # number of trees
                  shrinkage=0.005,             # shrinkage or learning rate,
                  # 0.001 to 0.1 usually work
                  interaction.depth=3,         # 1: additive model, 2: two-way interactions, etc.
                  bag.fraction = 0.5,          # subsampling fraction, 0.5 is probably best
                  train.fraction = 0.5,        # fraction of data for training,
                  # first train.fraction*N used for training
                  n.minobsinnode = 10,         # minimum total weight needed in each node
                  cv.folds = cv.fold,                # do 5-fold cross-validation
                  keep.data=TRUE,              # keep a copy of the dataset with the object
                  verbose=F)                # print out progress
      best.iter <- gbm.perf(gbm1,method="OOB", plot.it = F)
      data.predict = predict(gbm1, newdata = as.data.frame(x[-good.data,-j]), n.trees = best.iter)
      x[-good.data,j] = data.predict
      x[,j]
    })
  }
  
  return ( list (
    x=x,
    missing.matrix=missing.matrix
  ))
}



## ===============================================================
a <- gbmImpute(data2, max.iters = 2, cv.fold = 5, verbose=T)
names(a)
dim(a$x)
dim(data2)
sum(is.na(a$x))
