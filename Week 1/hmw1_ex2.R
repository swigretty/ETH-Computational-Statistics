#part a)
slopes <- c() #create an empty vector to store the slopes in
#create predictors
x   <- seq(1,40,1)  
#repeat 100 times
for (i in 1:100) {
  #simulate dependent variable
  y   <- 2*x+1+5*rnorm(length(x))
  reg <- lm(y~x) #fit a linear model
  slope  <- summary(reg)$coefficients[2]  #access the slope coefficient
  slopes <- c(slopes, slope)  #append it to vector
}

#part b)
par(mfrow=c(1,2)) #for two horisontal plots
hist(slopes, freq=FALSE, main = paste("Normaly distributed noise"))  #plot histogram of slopes
#append a column of ones to the design matrix, corresponding to the intercept
ones <- rep(1,length(x))
x <- cbind(ones, x)
#compute inverse of trnsopse(x)*x
inverse <- solve(t(x)%*%x)
#plot true density
lines(seq(1.8,2.3,by=0.01),dnorm(seq(1.8,2.3,by=0.01),mean=2,sd=5*sqrt(inverse[2,2])))

#part c)
#mean and standard error
mean <- mean(slopes)
std <- sd(slopes)

#part d)
x   <- seq(1,40,1)   #create predictors
#repeat 100 times
for (i in 1:100) {
  y   <- 2*x+1+5*(1-rchisq(length(x), df=1))/sqrt(2)
  reg <- lm(y~x)
  slope  <- summary(reg)$coefficients[2]  #access the slope
  slopes <- c(slopes, slope)  #append it to vector
}
hist(slopes, freq=FALSE, main = paste("Chi-squared noise"))  #plot histogram of slopes
lines(seq(1.8,2.3,by=0.01),dnorm(seq(1.8,2.3,by=0.01),mean=2,sd=5*sqrt(inverse[2,2])))
