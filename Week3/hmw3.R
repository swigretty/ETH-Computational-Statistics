##define the function
fun <- function(x) {
  return(x + 4*cos(7*x))
}
##set prediction points
x <- seq(-1, 1, length = 101)
m <- fun(x)
##set seed for reproducibility
set.seed(79)
## nw = Nadaraya-Watson, lp = Local Polynomial, ss = Smoothing Splines
nrep <- 1000
estnw <- estlp <- estss <- data <- matrix(0, nrow = 101, ncol = nrep)
##run simulations 
for(i in 1:nrep){
  ## Simulate y-values and store them
  y <- m + rnorm(length(x))
  data[,i] <- y
  ## Get estimates for the mean function m(x) in the points given by the vector x
  estnw[,i] <- ksmooth(x, y, kernel = "normal", bandwidth = 0.2, x.points = x)$y
  estlp[,i] <- predict(loess(y ~ x, span = 0.2971339), newdata = x)
  estss[,i] <- predict(smooth.spline(x,y,spar = 0.623396), x = x)$y
}
##estimate bias for the three estimators
meansnw <- apply(estnw, 1, mean)
biasnw <- meansnw - m
meanslp <- apply(estlp, 1, mean)
biaslp <- meanslp - m
meanss <- apply(estss, 1, mean)
biasss <- meanss - m
biases <- cbind(biasnw, biaslp, biasss)

##estimate variance for the three estimators
varnw <- apply(estnw, 1, var)
varlp <- apply(estlp, 1, var)
varss <- apply(estss, 1, var)
variances <- cbind(varnw, varlp, varss)

##estimate MSE for the three estimators
msenw <- apply((estnw - m)^2, 1, mean)
mselp <- apply((estlp - m)^2, 1, mean)
msess <- apply((estss - m)^2, 1, mean)
mse <- cbind(msenw, mselp, msess)

##plot everything
par(mfrow=c(1,3)) #for three horisontal plots
##biases
matplot(x,biases, pch = 1:3, col = 1, main = "Biases")
legend("top", inset=.05, legend=c("nw", "lp", "ss"), pch=1:3)
##variances
matplot(x,variances, pch = 1:3, col = 1, main = "Variances")
legend("top", inset=.05, legend=c("nw", "lp", "ss"), pch=1:3)
##MSEs
matplot(x,mse, pch = 1:3, col = 1, main = "MSEs")
legend("top", inset=.05, legend=c("nw", "lp", "ss"), pch=1:3)

## Plot MSE decomposition chcek
par(mfrow=c(1,1))
matplot(x,(mse - biases^2 - variances), pch = 1:3, col = 1,
          main = "MSE decomposition chcek", ylab = "MSE - bias^2 - var")
legend("bottom", inset=.05, legend=c("nw", "lp", "ss"), pch=1:3, horiz = TRUE)

##part b)
##calculate smoothing matrices for all three estimators
## nw = Nadaraya-Watson, lp = Local Polynomial, ss = Smoothing Splines
Snw <- Slp <- Sss <- matrix(0, nrow = 101, ncol = 101)
In <- diag(101) ## identity matrix
for(j in 1:101){
  y <- In[,j]
  Snw[,j] <- ksmooth(x, y, kernel = "normal", bandwidth = 0.2, x.points = x)$y
  Slp[,j] <- predict(loess(y ~ x, span = 0.2971339), newdata = x)
  Sss[,j] <- predict(smooth.spline(x,y,spar = 0.623396), x = x)$y
}
##calculate degrees of freedom in the estimators
dfnw <- sum(diag(Snw))
dflp <- sum(diag(Slp))
dfss <- sum(diag(Sss))
##calculate standard error estimates while running simulations
senw <- selp <- sess <- matrix(0, nrow = 101, ncol = nrep)
##run simulations and count the coverage rates of confidence intervals
coveragenw <-coveragelp <- coveragess <- 0 
coverageFullnw <-coverageFulllp <- coverageFullss <- 0 
for(i in 1:nrep){
  ## Get estimate for standard error
  sigma2nw <- sum((estnw[,i] - data[,i])^2) / (length(data[,i]) - dfnw)
  sigma2lp <- sum((estlp[,i] - data[,i])^2) / (length(data[,i]) - dflp)
  sigma2ss <- sum((estss[,i] - data[,i])^2) / (length(data[,i]) - dfss)
  senw[,i] <- sqrt(sigma2nw * diag(Snw))
  selp[,i] <- sqrt(sigma2lp * diag(Slp))
  sess[,i] <- sqrt(sigma2ss * diag(Sss))
  ## Check if estimate at 0.5 belongs to 95% confidence interval
  ## Note x[76] = 0.5. Calculate boundaries of bias adjusted confidence interval
  leftnw <- estnw[76,i] - 1.96*senw[76,i] - biasnw[76]
  rightnw <- estnw[76,i] + 1.96*senw[76,i] - biasnw[76]
  leftlp <- estlp[76,i] - 1.96*selp[76,i] - biaslp[76]
  rightlp <- estlp[76,i] + 1.96*selp[76,i] - biaslp[76]
  leftss <- estss[76,i] - 1.96*sess[76,i] - biasss[76]
  rightss <- estss[76,i] + 1.96*sess[76,i] - biasss[76]
  ##nw
  if ((m[76] >= leftnw) & (m[76] <= rightnw)) {
    coveragenw = coveragenw + 1
  }
  ##lp
  if ((m[76] >= leftlp) & (m[76] <= rightlp)) {
    coveragelp = coveragelp + 1
  }
  ##ss
  if ((m[76] >= leftss) & (m[76] <= rightss)) {
    coveragess = coveragess + 1
  }
  ## Check if estimate at all points belongs to 95% confidence interval
  bandcovernw <- bandcoverlp <- bandcoverss <- TRUE
  for (j in 1:101) {
    leftnw <- estnw[j,i] - 1.96*senw[j,i] - biasnw[j]
    rightnw <- estnw[j,i] + 1.96*senw[j,i] - biasnw[j]
    leftlp <- estlp[j,i] - 1.96*selp[j,i] - biaslp[j]
    rightlp <- estlp[j,i] + 1.96*selp[j,i] - biaslp[j]
    leftss <- estss[j,i] - 1.96*sess[j,i] - biasss[j]
    rightss <- estss[j,i] + 1.96*sess[j,i] - biasss[j]
    ##nw
    if ((m[j] < leftnw) || (m[j] > rightnw)) {
      bandcovernw <- FALSE
    }
    ##lp
    if ((m[j] < leftlp) || (m[j] > rightlp)) {
      bandcoverlp <- FALSE
    }
    ##ss
    if ((m[j] < leftss) || (m[j] > rightss)) {
      bandcoverss <- FALSE
    }
  }
  ## If all points in band increment counters
  if (bandcovernw) {
    coverageFullnw = coverageFullnw + 1
  }
  if (bandcoverlp) {
    coverageFulllp = coverageFulllp + 1
  }
  if (bandcoverss) {
    coverageFullss = coverageFullss + 1
  }
}

## For part c)
## copy-pasted the same code as above with changed x
##set prediction points
x <- sort(c(0.5, -1 + rbeta(50, 2, 2), rbeta(50, 2, 2)))
m <- fun(x)
##set seed for reproducibility
set.seed(79)
## nw = Nadaraya-Watson, lp = Local Polynomial, ss = Smoothing Splines
nrep <- 1000
estnw <- estlp <- estss <- data <- matrix(0, nrow = 101, ncol = nrep)
##run simulations 
for(i in 1:nrep){
  ## Simulate y-values and store them
  y <- m + rnorm(length(x))
  data[,i] <- y
  ## Get estimates for the mean function m(x) in the points given by the vector x
  estnw[,i] <- ksmooth(x, y, kernel = "normal", bandwidth = 0.2, x.points = x)$y
  estlp[,i] <- predict(loess(y ~ x, span = 0.37614), newdata = x)
  estss[,i] <- predict(smooth.spline(x,y,spar = 0.79424), x = x)$y
}
##estimate bias for the three estimators
meansnw <- apply(estnw, 1, mean)
biasnw <- meansnw - m
meanslp <- apply(estlp, 1, mean)
biaslp <- meanslp - m
meanss <- apply(estss, 1, mean)
biasss <- meanss - m
biases <- cbind(biasnw, biaslp, biasss)

##estimate variance for the three estimators
varnw <- apply(estnw, 1, var)
varlp <- apply(estlp, 1, var)
varss <- apply(estss, 1, var)
variances <- cbind(varnw, varlp, varss)

##estimate MSE for the three estimators
msenw <- apply((estnw - m)^2, 1, mean)
mselp <- apply((estlp - m)^2, 1, mean)
msess <- apply((estss - m)^2, 1, mean)
mse <- cbind(msenw, mselp, msess)

##plot everything
par(mfrow=c(1,3)) #for three horisontal plots
##biases
matplot(x,biases, pch = 1:3, col = 1, main = "Biases")
legend("top", inset=.05, legend=c("nw", "lp", "ss"), pch=1:3)
##variances
matplot(x,variances, pch = 1:3, col = 1, main = "Variances")
legend("top", inset=.05, legend=c("nw", "lp", "ss"), pch=1:3)
##MSEs
matplot(x,mse, pch = 1:3, col = 1, main = "MSEs")
legend("top", inset=.05, legend=c("nw", "lp", "ss"), pch=1:3)

## Plot MSE decomposition chcek
par(mfrow=c(1,1))
matplot(x,(mse - biases^2 - variances), pch = 1:3, col = 1,
        main = "MSE decomposition chcek", ylab = "MSE - bias^2 - var")
legend("bottom", inset=.05, legend=c("nw", "lp", "ss"), pch=1:3, horiz = TRUE)

##calculate smoothing matrices for all three estimators
## nw = Nadaraya-Watson, lp = Local Polynomial, ss = Smoothing Splines
Snw <- Slp <- Sss <- matrix(0, nrow = 101, ncol = 101)
In <- diag(101) ## identity matrix
for(j in 1:101){
  y <- In[,j]
  Snw[,j] <- ksmooth(x, y, kernel = "normal", bandwidth = 0.2, x.points = x)$y
  Slp[,j] <- predict(loess(y ~ x, span = 0.2971339), newdata = x)
  Sss[,j] <- predict(smooth.spline(x,y,spar = 0.623396), x = x)$y
}
##calculate degrees of freedom in the estimators
dfnw <- sum(diag(Snw))
dflp <- sum(diag(Slp))
dfss <- sum(diag(Sss))
##calculate standard error estimates while running simulations
senw <- selp <- sess <- matrix(0, nrow = 101, ncol = nrep)
##run simulations and count the coverage rates of confidence intervals
coveragenw <-coveragelp <- coveragess <- 0 
coverageFullnw <-coverageFulllp <- coverageFullss <- 0 
for(i in 1:nrep){
  ## Get estimate for standard error
  sigma2nw <- sum((estnw[,i] - data[,i])^2) / (length(data[,i]) - dfnw)
  sigma2lp <- sum((estlp[,i] - data[,i])^2) / (length(data[,i]) - dflp)
  sigma2ss <- sum((estss[,i] - data[,i])^2) / (length(data[,i]) - dfss)
  senw[,i] <- sqrt(sigma2nw * diag(Snw))
  selp[,i] <- sqrt(sigma2lp * diag(Slp))
  sess[,i] <- sqrt(sigma2ss * diag(Sss))
  ## Check if estimate at 0.5 belongs to 95% confidence interval
  ## Note x[76] = 0.5. Calculate boundaries of bias adjusted confidence interval
  leftnw <- estnw[76,i] - 1.96*senw[76,i] - biasnw[76]
  rightnw <- estnw[76,i] + 1.96*senw[76,i] - biasnw[76]
  leftlp <- estlp[76,i] - 1.96*selp[76,i] - biaslp[76]
  rightlp <- estlp[76,i] + 1.96*selp[76,i] - biaslp[76]
  leftss <- estss[76,i] - 1.96*sess[76,i] - biasss[76]
  rightss <- estss[76,i] + 1.96*sess[76,i] - biasss[76]
  ##nw
  if ((m[76] >= leftnw) & (m[76] <= rightnw)) {
    coveragenw = coveragenw + 1
  }
  ##lp
  if ((m[76] >= leftlp) & (m[76] <= rightlp)) {
    coveragelp = coveragelp + 1
  }
  ##ss
  if ((m[76] >= leftss) & (m[76] <= rightss)) {
    coveragess = coveragess + 1
  }
  ## Check if estimate at all points belongs to 95% confidence interval
  bandcovernw <- bandcoverlp <- bandcoverss <- TRUE
  for (j in 1:101) {
    leftnw <- estnw[j,i] - 1.96*senw[j,i] - biasnw[j]
    rightnw <- estnw[j,i] + 1.96*senw[j,i] - biasnw[j]
    leftlp <- estlp[j,i] - 1.96*selp[j,i] - biaslp[j]
    rightlp <- estlp[j,i] + 1.96*selp[j,i] - biaslp[j]
    leftss <- estss[j,i] - 1.96*sess[j,i] - biasss[j]
    rightss <- estss[j,i] + 1.96*sess[j,i] - biasss[j]
    ##nw
    if ((m[j] < leftnw) || (m[j] > rightnw)) {
      bandcovernw <- FALSE
    }
    ##lp
    if ((m[j] < leftlp) || (m[j] > rightlp)) {
      bandcoverlp <- FALSE
    }
    ##ss
    if ((m[j] < leftss) || (m[j] > rightss)) {
      bandcoverss <- FALSE
    }
  }
  ## If all points in band increment counters
  if (bandcovernw) {
    coverageFullnw = coverageFullnw + 1
  }
  if (bandcoverlp) {
    coverageFulllp = coverageFulllp + 1
  }
  if (bandcoverss) {
    coverageFullss = coverageFullss + 1
  }
}