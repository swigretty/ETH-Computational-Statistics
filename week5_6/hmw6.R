###################################################
### TASK a)
###################################################

### Function for fitting coefficients to data
lmcoefs <- function(data, ind) {
  d <- as.matrix(data)[ind,,drop=FALSE]
  coef(lm.fit(cbind(1, d[,c("x2","x3")]), d[,"y"]))
}

obst.est <- function(data, B)
{
  len <- nrow(data)
  # Create matrix of indexes of dimension (B x len) where
  # each colum corresponds to a bootstrap sample
  ind <- replicate(B, sample(1:len, len, replace = TRUE))
  # Bootstrap estimations of regression coefficients using
  # the B bootstrap samples
  apply(ind, 2, lmcoefs, data = data)
}

obst.ci <- function(bst.pars, data, alpha)
{
  ## Estimate regression parameters without using bootstrap
  reg.pars <- lmcoefs(data, 1:nrow(data))
  ## Calculate empirical quantiles of the bootstrap distribution
  qt <- apply(bst.pars, 1, quantile, probs = c(1 - alpha/2, alpha/2), names = FALSE)
  ## Return vector of bootstrap confidence intervals
  ## (cf. Formula (5.5) in the lecture notes)
  2*reg.pars - t(qt)
}

###################################################
### TASK b)
###################################################
set.seed(84)
num_simulations <- 1000
B <- 1000
x2 <- rep(1:5, 5)
x3 <- rep(1:5, each = 5)
grid <- cbind(x2,x3)
b <- cbind(1,-2,3)
alpha <- 0.05

### Normal
classical_coverage <- c(0,0,0)
bst_coverage <- c(0,0,0)
for (j in 1:num_simulations) {
  print(j)
  y <- rep(1,nrow(grid)) - 2*grid[,1] + 3*grid[,2] + rnorm(nrow(grid))
  data <- cbind(y, grid)
  # Construct classical confidence interval and check coverage
  fit <- lm(y ~ x2 + x3, data = data.frame(data))
  classical <- confint(fit, level = 1 - alpha)
  for (i in 1:3) {
    if ((classical[i,1] <= b[i]) && (classical[i,2] >= b[i])) {
      classical_coverage[i] <- classical_coverage[i] + 1
    } 
  }
  
  # Construct bootstrap confidence interval and check coverage
  bst.est <- obst.est(data,B)
  bst.int <- obst.ci(bst.est, data, alpha)
  for (i in 1:3) {
    if ((bst.int[i,1] <= b[i]) && (bst.int[i,2] >= b[i])) {
      bst_coverage[i] <- bst_coverage[i] + 1
    } 
  }
}

