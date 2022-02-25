## part a)
## read data
bmwlr <- scan("http://stat.ethz.ch/Teaching/Datasets/bmw.dat")
y <- bmwlr^2

## plot autocorolation for Xt and Xt^2
par(mfrow=c(1,2)) #for two vertical plots
acf(bmwlr, main = "Autocorrelation of X(t)")
acf(y, main = "Autocorrelation of X(t)^2")

## part b)
x <- sort(bmwlr)
y <- x^2
x <- x[1:999]
y <- y[2:1000]
par(mfrow=c(1,1))
plot(x,y, main = "Observed values under model", ylab="Y(t)", xlab="X(t-1)^2")
## fit local polynomial to data and chcek degrees of freedom
bmwloess <- loess(y ~ x)
df <- bmwloess$trace.hat
## fit a smoothing spline with found degrees of freedom
bmwsmooth <- smooth.spline(x,y,df=df)
## Function delta.dgf returns the degrees of freedom for the Nadaraya-Watson kernel
## regression with bandwidth 'h' minus the degrees of freedom 'dgf' to be matched
delta.dgf <- function (h, x, dgf){
  hatMat(x, trace = TRUE,
         pred.sm = function(x,y,...) ksmooth(sort(x), y, "normal", x.points=x, ...)$y,
         bandwidth = h) - dgf
}
library(sfsmisc)
bandwidth.nw <- uniroot(delta.dgf, c(3,4), x = x, dgf = df, tol = 0.01)$root
## fit a smoothing spline with found degrees of freedom
bmwnw <- ksmooth(x, y, kernel = "normal", bandwidth = bandwidth.nw, x.points = x)

## create residual plots
par(mfrow=c(1,3)) #for three vertical plots
plot(x, y - bmwnw$y, main = "Residual plot for NW", ylab = "residuals")
plot(x, y - fitted(bmwloess), main = "Residual plot for LP", ylab = "residuals")
plot(x, y - fitted(bmwsmooth), main = "Residual plot for SS", ylab = "residuals")
## create estimate plots
par(mfrow=c(1,3)) 
plot(x,bmwnw$y, main = "Estimate plot for NW", ylab = "estimates")
plot(x,fitted(bmwloess), main = "Estimate plot for LP", ylab = "estimates")
plot(x,fitted(bmwsmooth), main = "Estimate plot for SS", ylab = "estimates")

## create Tukey-Anscombe plots
par(mfrow=c(1,3)) #for three vertical plots
plot(bmwnw$y, y - bmwnw$y, main = "TA plot for NW", ylab = "residuals", xlab = "Fitted values")
plot(fitted(bmwloess), y - fitted(bmwloess), main = "TA plot for LP", ylab = "residuals", xlab = "Fitted values")
plot(fitted(bmwsmooth), y - fitted(bmwsmooth), main = "TA plot for SS", ylab = "residuals", xlab = "Fitted values")

## part c)
## fit data with lokerns and glkerns
library(lokern)
local <- lokerns(x, y, x.out=x)
global <- glkerns(x, y, x.out=x)
par(mfrow=c(1,2)) #for three vertical plots
plot(x, local$est, main = "Lokern estimate", ylab = "estimate")
plot(x, global$est,main = "Glkern estimate", ylab = "estimate")
par(mfrow=c(1,2)) #for three vertical plots
plot(x, (y-local$est),main = "Lokern residuals", ylab = "residuals")
plot(x, (y-global$est),main = "Glkern residuals", ylab = "residuals")
par(mfrow=c(1,1))
plot(x, (local$bandwidth), ylab = "Bandwidth", 
     main = "Local and global bandwidths", pch = 1)
abline(h = (global$bandwidth), lty = 4)
legend("topleft", inset=.1, legend=c("global"), lty = 4)
legend("top", inset=.1, legend=c("local"), pch = 1)
rug(x)
