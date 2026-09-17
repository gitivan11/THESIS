# Produce figure with
#   column 1: population and samples
#   column 2: sampling distribution for xbar, & bootstrap distributions
#   column 3: sampling distribution for t,    & bootstrap distributions
# for the mean from exponential populations, to show acceleration.

library(resample)
source("functions.R")

#--------------------------------------------------

# Work with exponential distributions
mu <- 1

#--------------------------------------------------

myColors <- c(4, 2, 3, 5, 6) # blue first, then red green cyan magenta

PlotSample <- function(data, i, lineMu = TRUE, star = FALSE, lineM = NULL,
                       col = myColors[1]){
  # lineMu: logical, add a vertical line at mu
  # star:  logical, label with xbar* instead of xbar
  # lineM: NULL or numeric; add a vertical line there

  hist(data, breaks = Breaks, col = col,
       probability = TRUE, axes = FALSE, new = FALSE, border = 0,
       xlab = "", xlim = Xlim, ylim = Ylim, xaxs = "i", yaxs = "i", main = "")
  axis(side = 1, at = Xlim)
  axis(side = 1, at = mean(data),
       if(star) expression(bar(x)^"*") else expression(bar(x)))
  if(lineMu) {
    segments(mu, y0 = 0, y1 = segmentTop)
    axis(side = 1, at = mu, "")
  }
  if(length(lineM)) {
    segments(lineM, y0 = 0, y1 = segmentTop, col = "red", lty = 2)
    axis(side = 1, at = lineM, "", col = "red")
  }
  axis(1, mu, expression(mu))
  text(Xlim[2], Ylim[2]*.8, adj = 1, paste("Sample", i))
  invisible(NULL)
}


#--------------------------------------------------
# Preliminary investigation: is there higher correlation between mean & sd
# in small or large samples?
f <- function(n, reps = 10000, returnValues = FALSE, returnData = FALSE) {
  x <- matrix(rexp(n*reps), n)
  if(returnData)
    return(x)
  temp <- list(mean = colMeans(x), sd = colStdevs(x))
  plot(temp$mean, temp$sd, pch = if(reps > 1000) "." else 1)
  if(returnValues)
    return(data.frame(temp, max = apply(x, 2, max)))
  c(meanOfMean = mean(temp$mean),
    sdOfMean = sd(temp$mean),
    meanOfSd = mean(temp$sd),
    sdOfSd = sd(temp$sd),
    cor = cor(temp$mean, temp$sd))
}
set.seed(0)
f(10)
f(100)
f(1000)
f(10000)
# correlations .757, .719, .713, .704
# So mild decrease with sample size
f(50)

# I'll use an exponential distribution (truncate picture at 6, don't use
# samples with larger values), n = 50. Select datasets with average sd for
# their mean.
n <- 50
set.seed(0)
Asummary <- f(n, returnValues = FALSE)
set.seed(0)
Astats <- f(n, returnValues = TRUE)
temp <- lm(sd ~ mean, data = Astats)
abline(temp)
Astats$resid <- resid(temp)
i <- order(Astats$mean)
Astats[head(i), ]
# obs 4744 has smallest mean and low resid
Astats[tail(i), ]
# obs 6091 has 5th largest mean, best resid, and largest value 5.55
K <- length(i)
which.min(abs(Astats$resid) + abs(Astats$mean - Astats$mean[i[.15*K]]))
which.min(abs(Astats$resid) + abs(Astats$mean - mu))
which.min(abs(Astats$resid) + abs(Astats$mean - Astats$mean[i[.85*K]]))
# use those to define Adata
set.seed(0)
Adata <- f(n, returnData = TRUE)[, c(4744, 9187, 1282, 199, 6091)]
points(colMeans(Adata), colStdevs(Adata), col = 2)

# Use Adata below. Turn it into a data frame, so can use sapply
Adata <- data.frame(Adata)
Adata.mean <- colMeans(Adata)
Adata.sd <- colStdevs(Adata)
B6 <- rep(c(10^3, 10^4), each = 3)
bootstraps <- lapply(Adata, bootstrap, statistic = mean, seed = 1, R = 10^4)
# bootstraps2 <- lapply(1:6,
#                       function(seed)
#                       bootstrap(data = Adata[, 4], statistic = mean,
#                                 seed = seed+1, R = B6[seed]))
# Would use bootstraps2 if I were doing more bootstrap distributions
# with different B. But I'll do t statistics instead.

tstat <- function(i) {
  # Run a bootstrap t
  data <- Adata[, i]
  meani <- mean(data)
  result <- bootstrap(data,
                      function(data) sqrt(n) * (mean(data) - meani) / sd(data),
                      R = 10^4, seed = i)
  result$observed <- NA
  result
}

bootstraps3 <- lapply(1:5, tstat)
hist(bootstraps3[[1]], xlim = c(-7.1, 3.4))
hist(bootstraps3[[2]], xlim = c(-7.1, 3.4))
hist(bootstraps3[[3]], xlim = c(-7.1, 3.4))
hist(bootstraps3[[4]], xlim = c(-7.1, 3.4))
hist(bootstraps3[[5]], xlim = c(-7.1, 3.4))
# xlim from -5 to 3.4 is fine. Can move mass lower than -5 to -5
bootstraps3b <- lapply(bootstraps3,
                      function(x) {
                        r <- x$replicates
                        r[r < -5] <- -5
                        x$replicates <- r
                        x})

hist(bootstraps3b[[1]], xlim = c(-5, 3.4))
hist(bootstraps3b[[2]], xlim = c(-5, 3.4))
hist(bootstraps3b[[3]], xlim = c(-5, 3.4))
hist(bootstraps3b[[4]], xlim = c(-5, 3.4))
hist(bootstraps3b[[5]], xlim = c(-5, 3.4))

#--------------------------------------------------
# Simulation trick to find the CDF of the t-statistic from exponential

# Variance Reduction Trick: xbar and (x / xbar) are independent.
# Draw samples of v = (x/xbar), then integrate over xbar
# xbar is gamma(shape = n, rate = n)
# Let V = S / xbar. Note that max(V) = sqrt(n)
# P(sqrt(n) (Xbar - mu)/S < a)
# = P(sqrt(n) (Xbar - mu)/(Xbar V) < a)  # monotone increasing in Xbar
# = P(Xbar - mu < Xbar V a / sqrt(n))
# = P(Xbar (1 - V a / sqrt(n)) < mu)
# = P(Xbar < mu / (1 - V a / sqrt(n))) # as long as (1-...) is positive
# = E_V( P(Xbar < mu / (1 - V a / sqrt(n)) | V))

xlim3 <- c(-5, 4) # xlim for column 3
x3 <- seq(from = xlim3[1], to = xlim3[2], length = 181)
set.seed(5)
reps <- 100 # increase this to 10K later
temp <- matrix(rexp(n*reps), n)
V <- colStdevs(temp) / colMeans(temp) # standardized sds
f <- function(V, x3) pgamma(mu / (1 - V * x3 / sqrt(n)), shape = n, rate = n)
tCDF <- colMeans(outer(V, x3, f))
plot(x3, tCDF, xlim = xlim3, type = "l", lwd = 2)
lines(x3, pt(x3, df = n-1))

tDensity <- diff(tCDF) / diff(x3)
x3b <- (x3[-1] + x3[-length(x3)])/2
plot(x3b, tDensity, xlim = xlim3, type = "l", lwd = 2)
lines(x3, dt(x3, df = n-1))

# Find .025 and .975 quantiles
tQuantiles <- approx(tCDF, x3, c(.025, .975))$y



#--------------------------------------------------

Xlim <- c(0, 6)
Breaks <- seq(Xlim[1], Xlim[2], length = 25) # width for each is 6/24 = .25
sapply(Adata, function(x) max(hist(x, breaks = Breaks, plot = FALSE)$counts))
# max is 18, for first
# density = fraction / barwidth = count / (n*.25)
Ylim <- c(0, 18 /(n*.25) + .01)
XlimMean <- range(Adata.mean) + c(-.5, .5)
YlimMean <- c(0, dnorm(0, sd = min(Adata.sd/sqrt(n))))
segmentTop <- Ylim * .9


#--------------------------------------------------
#pdf("../figures/distributionGraph4.pdf", height = 6, width = 6)

par(mfcol = c(6, 3), cex = .5, mar = c(4.1, .5, .1, .5))

### Column 1: Exponential density, and 5 samples

x <- seq(Xlim[1], Xlim[2], length = 61)
plot(x, dexp(x), type = "l", axes = FALSE, xaxs = "i", yaxs = "i",
     xlim = Xlim, ylim = Ylim, lwd = 2, xlab = "", ylab = "")
abline(h = 0, v = 0)
axis(side = 1, at = mu, expression(mu))
segments(mu, 0, y1 = segmentTop)
# for(i in c(1, 2, 4, 5)) # additional densities, exp matching data means
#   lines(x, dexp(x, rate = 1 / Adata.mean[i]), lwd = .5)
#axis(1, at = Adata.mean[-3], labels = rep("", 4))
text(Xlim[2], Ylim[2]*.8, adj = 1, paste("Population"))

# 5 samples
for(i in 1:5)
  PlotSample(Adata[, i], i, col = myColors[i])


### Column 2: Sampling Distribution, and 5 bootstrap distributions

x <- seq(XlimMean[1], XlimMean[2], length = 121)
plot(x, dgamma(x, shape = n, rate = n),
     type = "l", axes = FALSE, xaxs = "i", yaxs = "i",
     xlim = XlimMean, ylim = YlimMean, lwd = 2, xlab = "", ylab = "")
abline(h = 0, v = 0)
axis(side = 1, at = mu, expression(mu))
abline(v = mu)
# for(i in c(1, 2, 4, 5)) # additional densities, exp matching data means
#   lines(x, dgamma(x, shape = n, rate = n / Adata.mean[i]), lwd = .5)
# axis(1, at = Adata.mean[-3], labels = rep("", 4))
# Need axis call last, not abline, or arrows and text don't show up.
text(XlimMean[2], YlimMean[2]*.7, adj = c(1,0), paste("Sampling\nDistribution"))
text(XlimMean[2], YlimMean[2]*.7, adj = c(1,1.5),
     expression("for " * bar(x)))



# 5 bootstrap distributions
for(i in 1:5) {
  hist(bootstraps[[i]], col = myColors[i],
       xlim = XlimMean, ylim = YlimMean, xlab = "", ylab = "",
       axes = FALSE, yaxs = "i", legend = FALSE)
  abline(v = mu)
  axis(side = 1, at = Adata.mean[i], expression(bar(x)))
  axis(side = 1, at = Xlim)
  if(bootstraps[[i]]$observed > mu) {
    text(XlimMean[1], YlimMean[2]*.75, adj = 0,
         paste("Bootstrap\ndistribution\nfor sample", i))
  } else {
    text(XlimMean[2], YlimMean[2]*.75, adj = 1,
         paste("Bootstrap\ndistribution\nfor sample", i))
  }
}


### Column 3: bootstrap t distribution

# True density for t statistic
ylim3 <- c(0, dt(0, n-1)*1.04)
plot(x3b, tDensity, xlim = xlim3, type = "l", lwd = 2, yaxs = "i",
     axes = FALSE, ylim = ylim3, xlab = "", ylab = "")
# lines(x3, dt(x3, df = n-1))) # for reference
axis(side = 1, at = c(-4, -2, 0, 2, 4))
abline(v=0)
segments(tQuantiles, 0, y1 = par("usr")[4]*.5, lty = 3)
text(xlim3[2], ylim3[2]*.7, adj = c(1,0), paste("Sampling\nDistribution"))
text(xlim3[2], ylim3[2]*.7, adj = c(1,1.5),
     expression("for " * t))

for(i in 1:5) {
  hist(bootstraps3b[[i]], yaxs = "i", axes = FALSE, col = myColors[i],
       xlim = xlim3, ylim = c(0, dt(0, n-1)*1.04), xlab = "", ylab = "",
       legend = FALSE, showObserved = FALSE)
  axis(side = 1, at = c(-4, -2, 0, 2, 4))
  abline(v=0)
  temp <- quantile(bootstraps3b[[i]], c(.025, .975))
  segments(temp, 0, y1 = par("usr")[4]*.5, lty = 3)
  legend("topright", legend = "Quantiles",
         lty = 3, box.col = 0)
  text(-5, .2, adj = 0,
       paste("Bootstrap t\ndistribution\nfor sample", i))
}
#dev.off()
