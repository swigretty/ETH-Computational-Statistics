### Code skeleton for Series 5, task 2

## Read in dataset
diabetes <-
  read.table("http://stat.ethz.ch/Teaching/Datasets/diabetes2.dat",
             header = TRUE)
reg <- diabetes[, c("Age", "C.Peptide")]
names(reg) <-     c("x",   "y")

## Sort values
reg <- reg[sort.list(reg$x), ]

###################################################
### TASK a)
###################################################

### Utility function for LOO cross-validation:

##' Calculates the LOO CV score for given data and regression prediction function
##'
##' @param reg.data: regression data; data.frame with columns 'x', 'y'
##' @param reg.fcn:  regr.prediction function; arguments:
##'                    reg.x: regression x-values
##'                    reg.y: regression y-values
##'                    x:     x-value(s) of evaluation point(s)
##'                  value: prediction at point(s) x
##' @return LOOCV score
loocv <- function(reg.data, reg.fcn)
{
  ## Help function to calculate leave-one-out regression values
  loo.reg.value <- function(i, reg.data, reg.fcn)
    return(reg.fcn(reg.data[,'x'][-i], reg.data[,'y'][-i], reg.data[,'x'])[i])
  
  ## Calculate LOO regression values using the help function above
  n <- nrow(reg.data)
  loo.values <- sapply(1:n, loo.reg.value, reg.data, reg.fcn)
  
  ## Calculate and return MSE
  return(sum((reg.data[,'y']-loo.values)^2)/n)
}


### Regression prediction function for NW kernel:
reg.fcn.nw <- function(reg.x, reg.y, x)
  ksmooth(reg.x, reg.y, x.point = x, kernel = "normal", bandwidth = 3.5)$y

### Plot
est.nw <- reg.fcn.nw(reg[,'x'],reg[,'y'],reg[,'x'])
plot(reg[,'x'],reg[,'y'])
lines(reg[,'x'], est.nw)

### Calculation of LOO CV-score for NW kernel estimator:
(cv.nw <- loocv(reg, reg.fcn.nw))

### Hat matrix "S.nw":
n <- nrow(reg)
Id <- diag(n)
S.nw <- matrix(0, n, n)
for (j in 1:n)
  S.nw[, j] <- reg.fcn.nw(reg[,'x'], Id[,j], reg[,'x'])

### Degrees of freedom (cf. Formula (3.6) in the lecture notes:
(df.nw <- sum(diag(S.nw)))

###################################################
### TASK b)
###################################################

### Regression prediction function for local polynomial:
reg.fcn.lp <- function(reg.x, reg.y, x)
  predict(loess(reg.y ~ reg.x, enp.target = df.nw, surface = "direct"), newdata = x)

### Plot
est.lp <- reg.fcn.lp(reg[,'x'],reg[,'y'],reg[,'x'])
plot(reg[,'x'],reg[,'y'])
lines(reg[,'x'], est.lp)

### Calculation of LOO CV-score for local polynomial estimator:
(cv.lp <- loocv(reg, reg.fcn.lp))

###################################################
### TASK c) and d)
###################################################

est.ss <- smooth.spline(reg[,'x'],reg[,'y'],df = df.nw, cv = TRUE)
est.ss$cv.crit
### Regression prediction function for local polynomial:
reg.fcn.ss <- function(reg.x, reg.y, x)
  predict(smooth.spline(reg.x,reg.y, spar = est.ss$spar), x = x)$y

### Plot
y <- reg.fcn.ss(reg[,'x'],reg[,'y'],reg[,'x'])
plot(reg[,'x'],reg[,'y'])
lines(reg[,'x'], y)

### Calculation of LOO CV-score for smoothing spline estimator:
(cv.ss <- loocv(reg, reg.fcn.ss))

### Automatic selection of degrees of freedom
est.ss.auto <- smooth.spline(reg[,'x'],reg[,'y'], cv = TRUE)
est.ss.auto$cv.crit
plot(reg[,'x'],reg[,'y'])
lines(reg[,'x'], predict(est.ss.auto, x=reg[,'x'])$y)

###################################################
### TASK e)
###################################################

### Mean fit function
reg.fcn.const <- function(reg.x, reg.y, x)
  rep(mean(reg.y), length(x))

### Plot
est.const <- reg.fcn.const(reg[,'x'],reg[,'y'],reg[,'x'])
plot(reg[,'x'],reg[,'y'])
lines(reg[,'x'], est.const)

### Calculation of LOO CV-score for constant fit estimator:
(cv.const <- loocv(reg, reg.fcn.lp))
