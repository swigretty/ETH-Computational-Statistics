set.seed(1)
##initialize vector of bandwidths
band <- c(0.02,0.1,0.3,0.6,1,1.5,2)
qualities <- c()
##iterat through bandwidths
for (j in 1:length(band)) {
  quality <- numeric(200)
  ##simulate 200 times
  for (k in 1:200) {
    ##generate data from mixture of gaussians
    data <- numeric(100)
    for(i in 1:100){
      p <- runif(1, min = 0, max = 1)
      if (p < 0.2) {
        data[i] <- rnorm(1, mean = 0, sd = sqrt(0.01))
      }
      else {
        data[i] <- rnorm(1, mean = 2, sd = 1)
      }
    }
    ##fit density
    ke <- density(data, bw = band[j], n = 61, from = -1, to = 5)
    ## Compute the true density at the given datapoints
    dmix <- 0.2 * dnorm(ke$x[-1], mean = 0, sd = sqrt(0.01)) +
      0.8 * dnorm(ke$x[-1], mean = 2, sd = 1)
    ## Take the mean of the squared differences
    quality[k] <- mean((ke$y[-1] - dmix)^2)
  }
  qualities <- cbind(qualities, mean(quality))
}