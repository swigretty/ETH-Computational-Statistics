#part a) read and plot data
airline <- scan("http://stat.ethz.ch/Teaching/Datasets/airline.dat")
time  <- seq(1:length(airline))
par(mfrow=c(1,2))
plot(time, airline, type = "l", main = "Passenger count \n over time", xlab = "Months", ylab = "Number of passengers")
plot(time, log(airline), type = "l",  main = "Logarithmic Passenger \n count over time", xlab = "Months", ylab = "Log(Number of passengers)")

#create month vector indicators
f1 <- seq(1:length(airline))
f2 <- rep(c(1,rep(0,11)),12)
f3 <- rep(c(rep(0,1),1,rep(0,10)),12)
f4 <- rep(c(rep(0,2),1,rep(0,9)),12)
f5 <- rep(c(rep(0,3),1,rep(0,8)),12)
f6 <- rep(c(rep(0,4),1,rep(0,7)),12)
f7 <- rep(c(rep(0,5),1,rep(0,6)),12)
f8 <- rep(c(rep(0,6),1,rep(0,5)),12)
f9 <- rep(c(rep(0,7),1,rep(0,4)),12)
f10 <- rep(c(rep(0,8),1,rep(0,3)),12)
f11 <- rep(c(rep(0,9),1,rep(0,2)),12)
f12 <- rep(c(rep(0,10),1,rep(0,1)),12)
f13 <- rep(c(rep(0,11),1),12)

#fit linear model
reg <- lm(log(airline) ~ f1+f2+f3+f4+f5+f6+f7+f8+f9+f10+f11+f12+f13-1)
summary(reg)

#plot
par(mfrow=c(1,2))
plot(time, fitted(reg),type = "l", main = "Regression fitted \n logarithmic passenger \n count over time", xlab = "Months", ylab = "Fitted Log(Number of passengers)")
plot(time, resid(reg), main = "Residuals", xlab = "Months", ylab = "Residuals of Log(Number of passengers)")
