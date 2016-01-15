### Series 8, task 1
## Read in dataset
heart <- read.table("http://stat.ethz.ch/Teaching/Datasets/heart.dat",
                   header = TRUE)

## Negative log likelihood function
neg.ll <- function(beta, data) {
  ## first calculate log of number of combinations
  log.combinations <- log(choose(data[,'m'],data[,'N']))
  ## calculate g(beta,x)
  g <- beta[1] + beta[2]*data[,'age']
  ## calculate second and third term
  second <- data[,'N']*g
  third <- data[,'m']*log(1+exp(g))
  -sum(log.combinations + second - third)
}

## Plot the negative log likelihood
beta0.grid <- seq(-10, 10, length = 101)
beta1.grid <- seq(-10, 10, length = 101)
neg.ll.values <- matrix(0, nrow = 101, ncol = 101)
for (i in 1:101) {
  for (j in 1:101) {
    beta <- c(beta0.grid[i], beta1.grid[j])
    neg.ll.values[i,j] <- neg.ll(beta, heart)
  }
}
contour(beta0.grid, beta1.grid, neg.ll.values)


## Fit a beta
fit <- glm(cbind(N, m - N) ~ age, family = binomial, data = heart)

## Optimize
optim(c(0, 0), neg.ll, data = heart)

## get a probability estimate
new.age <- 1:100
probability.est <- predict(fit, newdata = data.frame(age = new.age), type = "response")
plot(1:100,probability.est)


### Code skeleton for Series 8, task 2

## Read in dataset
ozone <- read.table("http://stat.ethz.ch/Teaching/Datasets/ozone.dat",
                    header = TRUE)

###################################################
### TASK a)
###################################################

ozone$logupo3 <- log(ozone$upo3)
d.ozone <- subset(ozone, select=-upo3)
pairs(d.ozone, pch = ".",gap = 0.1)

## delete outlier
out <- which.max(d.ozone[,"wdsp"])
d.ozone.e <- d.ozone[-out,]



###################################################
### TASK b)
###################################################

## package for formula
require(sfsmisc)

## Linear models
## fit 1 (polynomial of degree 1)
form1 <- logupo3~.
fit1 <- lm(form1, d.ozone.e)



## fits of degree 2 to 5
form2 <- wrapFormula(logupo3~., d.ozone.e, wrapString="poly(*,degree=2)")
fit2 <- lm(form2, d.ozone.e)

form3 <- wrapFormula(logupo3~., d.ozone.e, wrapString="poly(*,degree=3)")
fit3 <- lm(form3, d.ozone.e)

form4 <- wrapFormula(logupo3~., d.ozone.e, wrapString="poly(*,degree=4)")
fit4 <- lm(form4, d.ozone.e)

form5 <- wrapFormula(logupo3~., d.ozone.e, wrapString="poly(*,degree=5)")
fit5 <- lm(form5, d.ozone.e)



## GAM
require(mgcv)
gamForm <- wrapFormula(logupo3~.,d.ozone.e, wrapString="s(*)")
g1 <- gam(gamForm, data = d.ozone.e)


###################################################
### TASK c)
###################################################

## plot the fits

source("ftp://stat.ethz.ch/Teaching/maechler/CompStat/plotGAM.R") # to get p.gam()
p.gam(g1)
###################################################
### TASK d)
###################################################


## Mallows Cp function

Cp <- function(object,sigma){
  res<-residuals(object)
  n <- length(res)
  p <- n-object$df.residual
  SSE <- sum(res^2)
  SSE/sigma^2-n+2*p
}

## set sigma (use estimated sigma from fit5)
sigma <- sd(fit5$residuals)

## Calculate Mallows's Cp statistic for all 6 models
Cp(fit1, sigma)
Cp(fit2, sigma)
Cp(fit3, sigma)
Cp(fit4, sigma)
Cp(fit5, sigma)
Cp(g1, sigma)