### Ch. 5 Bootstrap
### =================

## Running example using correlation coefficient
## ---------------------------------------------

set.seed(1)
x <- rnorm(50)
y <- 1+x + rnorm(50,0,1)

plot(x,y, main = "n=50 -- \"The Data\" ")

(th.n <- cor(x,y))
## 0.6339331

## In order to skip the "long" computations below, you can get the saved results:
saveFile <- "~/Vorl/comput-statist/R/cor_boot.rda"
saveURL <- "http://stat.ethz.ch/Teaching/maechler/CompStat/cor_boot.rda"
try(print(load(url(saveURL)))) # "res"  "res.true"  "r.bseqvar"
## if that failed:
if(file.exists(saveFile))
    load(saveFile)
##  --------------
## --> and you can skip to the  " jump " (text) below

B <- 100000
res <- numeric(B)
set.seed(22)
for(i in 1:B) {
    if(i %% 10 == 0) cat(".", if(i %% 200 == 0) sprintf(" %6d\n", i))
    ind <- sample.int(length(x), replace = TRUE)
    xb <- x[ind]
    yb <- y[ind]
    res[i] <- cor(xb,yb)
    ##        ----------
}

## The "true" distribution (also by simulation): use  'N ~= oo (Infinity)'
##               (slightly slower than bootstrap above - can you see why?)
N <- 100000
res.true <- numeric(N)
set.seed(22)
for(i in 1:N) {
    if(i %% 10 == 0) cat(".", if(i %% 200 == 0) sprintf(" %6d\n", i))
    x <- rnorm(50)
    y <- 1+x + rnorm(50,0,1)
    res.true[i] <- cor(x,y)
}


## <--- jump to here

th.lab <- function(th, label = expression(hat(theta)[n]), col = "skyblue3", lwd = 3,
                   tck = 0.2, padj = if(grepl("\\<hat\\>",format(label))) -.1 else -.4)
    axis(1, at = th, labels = label, tck=tck, padj=padj, lwd=lwd,
         col=col, col.axis=col)

require(sfsmisc)# mult.fig() etc

## Compare Bootstrap world with Real World:
pdf.do(file = "hist-2.pdf", height = 7, width = 7)

mult.fig(mfrow = c(2,1))
br <- seq(0,1, by= 1/128)
hist(res,     breaks = br, main = "Bootstrap, B=100'000")
th.lab(th.n)
hist(res.true,breaks = br, main = "true distribution, simulated 100'000 times")
th.lab(sqrt(1/2), expression(theta == sqrt(2)/2))

pdf.end()

mean(res) ## 0.6241057
var(res)  ## 0.008215957

mean(res.true)## 0.7033774  [ N -> oo  :  Cov(X,Y) = 1/sqrt(2) = 0.7071 ]
var(res.true) ## 0.005343749

quantile(res,c(0.025,0.975)) #-- see confint etc, below
##      2.5%     97.5%
## 0.4170626 0.7713382

## Compute Var(T*_1 .. T*_B)  for increasing B  (B = 10,20, .... 100'000) :
B. <- B / 10
r.bseqvar <- numeric(B.)
for(i in 1:B.) {
    if(i %% 10 == 0) cat(".", if(i %% 200 == 0) sprintf(" %6d\n", i))
    r.bseqvar[i] <- var(res[1:(10*i)])
}

if(FALSE) ## and if you are MM:
   ## Save the results for fast reloading
   save(B, res, res.true, r.bseqvar, file = saveFile)

require("sfsmisc") # for mult.fig(), pdf.do(), ..

pdf.do(file = "var-boot-conv.pdf", width = 6, height = 6)

mult.fig(mfrow = c(1,2), main = expression("Bootstrap variance  " *
				           var(T^"*"[1],...,T^"*"[B])))
plot(10*(1:100),r.bseqvar[1:100], type = "l",
     xlab = "B", ylab = "bootstrap variance", main = "B = 10 .. 1000")
abline(h = var(res),      col = 2)
abline(h = var(res.true), col = 4,lty = 2)

## show the "zoom region" of the next plot
px <- par("usr")[1:2]; bgcol <- rgb(.7,.8,.9, alpha=.3)
rect(px[1], 0.007, px[2], 0.01, col = bgcol, border = "gray70", lwd=2)

plot(10*(1:10000),r.bseqvar,type = "l", ylim = c(0.007,0.01), xaxt="n",# <- better x-labels
     xlab = "B", ylab = "bootstrap variance",  main = "B = 10 .. 100'000")
axis(1, at = (x0 <- c(10000, 100000)), labels = sprintf("%d", x0))
abline(h = var(res),      col = 2)

if(FALSE)# show the same shaded area as on the left [FALSE: it's ugly]
rect(px[1], 0.007, px[2], 0.01, col = bgcol, border = "gray70", lwd=.5)

pdf.end()


## Quantiles:

(q95 <- quantile(res, pr = c(0.025, 0.975)))
##      2.5%     97.5%
## 0.4170626 0.7713382
dens.true <- density(res.true)

pdf.do(file = "hist-CI.pdf", width = 11.6, height = 8.2)# A4 quer

par(mfrow = c(1,1))
hist(res, breaks = "FD", # <- better than default ("Sturges") algorithm
     main = "histogram of B=100'000 bootstrapped correlation estimators",
     prob = TRUE, ylim = c(0, max(dens.true$y)),
     col = "gray88", border = "gray60", xlab = "")
th.lab(th.n, padj=0)

axis(1, at = q95, tck =  0.3, col = 2, labels=FALSE)
mtext("2.5% and 97.5%  quantiles", 1, at = mean(q95), line = -1.1, col = 2)

### Bootstrap confidence intervals -- now do the confidence interval correctly

(q95. <- quantile(res - th.n, pr = c(0.025, 0.975)))
##       2.5%      97.5%
## -0.2168706  0.1374051

(ci95 <- th.n - q95.[2:1])
##     97.5%      2.5%
## 0.4965280 0.8508037

axis(1, at = ci95, tck = 0.1,
     col="seagreen", col.axis="seagreen", lwd = 2, labels=FALSE)
mtext("correct  95% - C.I.", 1, at = mean(ci95), line = -2.2,
      col = "seagreen")

lines(dens.true, col="seagreen")

pdf.end()
