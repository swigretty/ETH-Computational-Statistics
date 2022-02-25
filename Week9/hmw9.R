### Series 9, Task 1
### cv function
cv <- function(fitfn, formula = logupo3 ~ . , data = d.ozone.es, ..., trace = TRUE)
{
  modFrame <- model.frame(formula, data = data)
  nc <- nrow(data)
  ssr <- 0
  if(trace) cat("  j = ")
  for(j in 1:nc) {
    if(trace) cat(if(j %% (nc %/% 10) == 1) paste(j, "") else ".")
    ## Fit without 'j' :
    fit <- fitfn(formula=formula, data = data[-j ,], ...)
    ## Evaluate at 'j' :
    ssr <- ssr + (model.response(modFrame)[j] - predict(fit, modFrame[j,]))^2
  }
  if(trace) cat("\n")
  ssr
}

## Read in dataset
ozone <- read.table("http://stat.ethz.ch/Teaching/Datasets/ozone.dat",
                    header = TRUE)
###################################################
### TASK a)
###################################################
d.ozone <- subset(transform(ozone, logupo3 = log(upo3)), select = -upo3)
d.ozone.e <- d.ozone[-which.max(d.ozone[,"wdsp"]),]
d.ozone.es <- d.ozone.e
d.ozone.es[,-10] <- scale(d.ozone.e[,-10])
###################################################
### TASK b)
###################################################
form <- logupo3~.
earth.mod <- earth(form, data = d.ozone.es, degree = 2)
plotmo(earth.mod, degree2 = FALSE)
plotmo(earth.mod, degree1 = FALSE)
###################################################
### TASK c)
###################################################
ppr.mod <- ppr(form, data = d.ozone.es, nterms = 4)
par(mfrow=c(1,1))
plot(ppr.mod)
###################################################
### TASK d)
###################################################
smoothers <- c("supsmu", "spline", "gcvspline")
num.ridges <- c(3,4,5)
deg.freed <- c(5,6,7)
models <- c()
cv.score <- c()
### for each set of parameters calculate cv score
### for each smoother method
for (sm in smoothers) {
  ### for each number of ridge functions
  for (nr in num.ridges) {
    ### if sm is spline try different DFs
    if (sm == "spline") {
      for (df in deg.freed) {
        cat("Smoother method: ", sm, ", number of ridge terms = ", nr,
            ", euqivalent degrees of freedom = ", df, ".\n")
        model <- ppr(form, data = d.ozone.es, sm.method = sm, nterms = nr, df = df)
        models <- cbind(models, model)
        score <- cv(ppr, sm.method = sm, nterms = nr, df = df)
        cv.score <- cbind(cv.score, score)
        cat("Score = ", score, "\n")
      }
    }
    else {    
      cat("Smoother method: ", sm,", number of ridge terms = ", nr, ".\n")
      model <- ppr(form, data = d.ozone.es, sm.method = sm, nterms = nr)
      models <- cbind(models, model)
      score <- cv(ppr, sm.method = sm, nterms = nr)
      cv.score <- cbind(cv.score, score)
      cat("Score = ", score, "\n")
    }
  }
}
###################################################
### TASK e)
###################################################
deg <- c(1,2,3,4)
for (d in deg) {
  cat("Degree = ", d,".\n")
  score <- cv(earth, degree = d)
  cat("Score = ", score, ".\n")
}

require(mgcv)
gamForm <- wrapFormula(logupo3~.,d.ozone.e, wrapString="s(*)")
gam.cv <- cv(gam)
###################################################
### TASK f)
###################################################
best <- ppr(form, data = d.ozone.es, sm.method = "spline", nterms = 5, df = 6)
best.cv <- cv(ppr, sm.method = "spline", nterms = 5, df = 6)
plotmo(best)
