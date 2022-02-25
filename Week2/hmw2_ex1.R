## Reading the dataset
url <- "http://stat.ethz.ch/Teaching/Datasets/mortality.csv"
mortality <- read.csv(url,header = TRUE)
mortality <- mortality[,-1]
##Create pairs plot using the splom() function of the lattice package
library(lattice)
splom(~mortality,pscales=0,cex=0.5)

mortal.full <- lm(Mortality ~. , data=mortality)
plot(1:length(resid(mortal.full)), resid(mortal.full))
plot(fitted.values(mortal.full), resid(mortal.full))

mortality[4] <- log(mortality[4])
mortality[7] <- log(mortality[7])
mortality[10] <- log(mortality[10])
mortality[13] <- log(mortality[13])
mortality[14] <- log(mortality[14])
splom(~mortality,pscales=0,cex=0.5)

## Fit the full model
mortal.full <- lm(Mortality ~. , data=mortality)
plot(1:length(resid(mortal.full)), resid(mortal.full))
plot(1:length(fitted.values(mortal.full)), fitted.values(mortal.full))
## Fit the empty model. This is not very useful in itself, but is required
## as a starting model for stepwise forward variable selection
mortal.empty <- lm(Mortality ~ 1, data = mortality)
## Backward elimination, starting from the full model
mortal.bw    <- step(mortal.full, direction = "backward")
## Forward selection, starting from the empty model
mortal.fw    <- step(mortal.empty, direction = "forward",
                         scope = list(upper=mortal.full,lower=mortal.empty))
## Loading the package for all-subsets regression
library(leaps)
## All subsets model choice, compare to the stepwise methods
mortal.alls  <- regsubsets(Mortality ~. , data=mortality)
## Load function to produce a nice figure of C_p versus p
source("ftp://stat.ethz.ch/Teaching/maechler/CompStat/cp-plot.R")
p.regsubsets(mortal.alls)
