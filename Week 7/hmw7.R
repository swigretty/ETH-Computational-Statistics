### Code skeleton for Series 7, task 2

## Read in dataset, set seed, load package
Iris <- iris[,c("Petal.Length","Petal.Width","Species")]
grIris <- as.integer(Iris[,"Species"])
set.seed(16)
library(MASS)

## Read n
n <- nrow(Iris)

## Utility functiom for plotting boundaries
predplot <- function(object, x, gr = grIris, main = "", lines.only=FALSE,
                     len = 42, colcont = "black", ...)
{
  ##  gr : the true grouping/class vector
  stopifnot(length(gr) == nrow(x))
  xp <- seq(min(x[,1]), max(x[,1]), length=len)
  yp <- seq(min(x[,2]), max(x[,2]), length=len)
  grid <- expand.grid(xp, yp)
  colnames(grid) <- colnames(x)[-3]
  Z <- predict(object, grid, ...)
  zp <- as.numeric(Z$class)
  zp <- Z$post[,3] - pmax(Z$post[,2], Z$post[,1])
  if(!lines.only)
    plot(x[,1], x[,2], col =gr, pch = gr,
         main = main,xlab=colnames(x)[1],ylab=colnames(x)[2])
  contour(xp, yp, matrix(zp, len),
          add = TRUE, levels = 0, drawlabels = FALSE, col=colcont)
  zp <- Z$post[,1] - pmax(Z$post[,2], Z$post[,3])
  contour(xp, yp, matrix(zp, len),
          add = TRUE, levels = 0, drawlabels = FALSE, col=colcont)
}
## Bootstrap size
B <- 1000




###################################################
### TASK a)
###################################################

## Use function lda to fit data
class_lda <- lda(Species~., Iris)

## Use function predplot to plot the boundaries
predplot(class_lda, Iris, main="Classification with LDA")

## Use function qda to fit data
class_qda <- qda(Species~., Iris)

## Use function predplot to plot the boundaries
predplot(class_qda, Iris, main="Classification with QDA")


###################################################
### TASKS b)
###################################################

## Create a random index matrix with either functions sample or sample.int to generate bootstrap
index <- matrix(sample(n, n*B, replace = TRUE), nrow=n, ncol=B)

## Initialize the list for LDA nad QDA fits
fit_lda <- vector("list",B)
fit_qda <- vector("list",B)


## Use both methods on the bootstrap samples
for(i in 1:B) {
  ind <- index[,i]
  fit_lda[[i]] <- lda(Species~., Iris[ind,])
  fit_qda[[i]] <- qda(Species~., Iris[ind,])
}

## Initialize the mu_hat bootstrap estimates
mu_hat_1 <- mu_hat_2 <- mu_hat_3 <- matrix(0,ncol=B,nrow=2)

## Determine the mu_hat bootstrap estimates

for(i in 1:B){
  mu_hat_temp <- fit_lda[[i]]$means
  mu_hat_1[,i] <- mu_hat_temp[1,]
  mu_hat_2[,i] <- mu_hat_temp[2,]
  mu_hat_3[,i] <- mu_hat_temp[3,]
}

## Plot the boostrapped estimators
## get borders of plot
min_x <- min(Iris[,1])
min_y <- min(Iris[,2])
max_x <- max(Iris[,1])
max_y <- max(Iris[,2])
plot(mu_hat_1[1,],mu_hat_1[2,], xlim = c(min_x, max_x), ylim = c(min_y, max_y))
points(mu_hat_2[1,],mu_hat_2[2,], col = "red", pch = 1)
points(mu_hat_3[1,],mu_hat_3[2,], col = "green", pch = 1)

###################################################
### TASK c)
###################################################

## Plot the bootstrapped boundaries estimates with LDA
predplot(class_lda, Iris,
         main = "Bootstrapped boundaries estimates with LDA")
for(i in 1:B){
  fit <- fit_lda[[i]]
  predplot(fit, Iris, lines.only= TRUE, colcont=adjustcolor("gray", 0.25))
}



## Plot the bootstrapped boundaries estimates with QDA
predplot(class_qda, Iris,
         main= "Bootstrapped boundaries estimates with QDA")
for(i in 1:B){
  fit <- fit_qda[[i]]
  predplot(fit, Iris, lines.only= TRUE, colcont=adjustcolor("gray", 0.25))
}

###################################################
### TASK d)
###################################################


## Initialize the errors

error_lda <- rep(0,B)
error_qda <- rep(0,B)

## Use the predict function to calculate the error
## Read help on predict.lda. Remember that logical
## FALSE/TRUE are treated as 0/1.

for(i in 1:B){
  error_lda[i] <- mean(predict(fit_lda[[i]],Iris)$class != Iris[,3])
  error_qda[i] <- mean(predict(fit_qda[[i]],Iris)$class != Iris[,3])
}

## Print the error
cat("Generalized error for LDA:",format(mean(error_lda),digits=4))
cat("Generalized error for QDA:",format(mean(error_qda),digits=4))


## Plot the boxplot of the errors

boxplot(cbind(error_lda, error_qda), range = 1.5)

###################################################
### TASK e)
###################################################


## Initialize the errors

error_lda_oob <- rep(0,B)
error_qda_oob <- rep(0,B)

## Use the predict function to calculate the error
## Read help on predict.lda. Remember that logical
## FALSE/TRUE are treated as 0/1.

for(i in 1:B){
  ind <- index[,i]
  error_lda_oob[i] <- mean(predict(fit_lda[[i]],Iris[-ind,])$class != Iris[-ind,3])
  error_qda_oob[i] <- mean(predict(fit_qda[[i]],Iris[-ind,])$class != Iris[-ind,3])
}

## Print the error
cat("Generalized error for LDA:",format(mean(error_lda_oob),digits=4))
cat("Generalized error for QDA:",format(mean(error_qda_oob),digits=4))


## Plot the boxplot of the errors

boxplot(cbind(error_lda_oob, error_qda_oob), range = 1.5)

