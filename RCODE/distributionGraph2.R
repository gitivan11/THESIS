# Produce figure with
#   column 1: population and samples
#   column 2: sampling distribution and bootstrap distributions
#   column 3: additional bootstrap distributions
# for the mean, small odd n.

library(resample)
source("functions.R")

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
  }
  if(length(lineM)) {
    segments(lineM, y0 = 0, y1 = segmentTop, col = "red", lty = 2)
    axis(side = 1, at = lineM, "", col = "red")
  }
  # axis(1, mu, expression(mu))
  text(Xlim[2], Ylim[2]*.8, adj = 1, paste("Sample", i))
  invisible(NULL)
}


#--------------------------------------------------
# Define population, limits, samples

#### Density for Population
# Population is N(0, 2^2)

n <- 9
mu <- 0
sd1 <- 2

# Some values for use in plotting
Xlim <- c(-6, 6)
Breaks <- seq(Xlim[1], Xlim[2], length = 21)
maxBins <- 3 # Only use samples with short histograms to avoid clipping

# density = fraction / barwidth = count / (n*14/20), max = maxBins/(...)
Ylim <- c(0, (maxBins / n) / (diff(Xlim)/(length(Breaks)-1)))
segmentTop <- .9 * Ylim[2]

Ylim*sqrt(n) # This is taller than needed
YlimMean <- 1.2 * c(0, dnorm(0, sd = sd1/sqrt(n)))
XlimMean <- c(-3, 3)




# Generate some random samples; I'll pick some with varying mu and sigma
set.seed(0)
A <- lapply(1:300, function(i, n) rnorm(n, sd = sd1), n = n)
A1means <- sapply(A, mean)
A1sd <- sapply(A, sd)
A1skew <- sapply(A, function(x) mean(x-mean(x))^3) / A1sd^3

if(FALSE) {
  plot(A1means, A1sd)
  abline(h = sd1)
  # Don't want any with values outside Xlim
  bad <- (sapply(A, min) < Xlim[1]) | (sapply(A, max) > Xlim[2])
  PointsSubset(A1means, A1sd, bad, col = "yellow")
  # Don't want any samples with a too-tall bin
  A[bad] <- lapply(A[bad], function(x) mu)
  temp <- lapply(A, hist, breaks = Breaks, plot = FALSE)
  tallBin <- sapply(temp, function(x) max(x$counts))
  table(tallBin) # table for number in tallest bin
  PointsSubset(A1means, A1sd, tallBin > maxBins, col = "yellow")
  # Show some with large and small skewness
  PointsSubset(A1means, A1sd, A1skew < quantile(A1skew, .1), pch = ".")
  PointsSubset(A1means, A1sd, A1skew > quantile(A1skew, .9), pch = "+")
  text(A1means[ii], A1sd[ii], ii, col = "green") # show current selection
  identify(A1means, A1sd, n = 5)
  # Save interesting ones as ii. Want first to have average sd & moderately
  # extreme mean.
}
# ii <- as.integer(c(157, 58, 86, 237, 288))
ii <- as.integer(c(157, 58, 68, 295, 288))
first <- ii[1] # low mean, avg sd

B6 <- rep(c(10^3, 10^4), each = 3)
bootstraps2 <- lapply(ii,
                     function(i) bootstrap(A[[i]], mean, seed = i, R = 10^4))
bootstraps3 <- lapply(1:6,
                     function(i) bootstrap(A[[first]], mean, seed = i, R = B6[i]))


#--------------------------------------------------
# Make distributionGraph2

# Either 5 or 4 datasets (below the population)
nDatasets <- 4  # or 5

#pdf("../figures/distributionGraph2.pdf", height = 1+nDatasets, width = 6)
par(mfcol = c(1+nDatasets, 3), cex = .5, mar = c(4.1, .5, .1, .5))

### Column 1: population density, then 5 samples
x <- seq(from = Xlim[1], to = Xlim[2], length = 151)
plot(x, dnorm(x, sd = sd1), xlab = "", ylab = "",
     type = "l", ylim = Ylim, xlim = Xlim,
     axes = FALSE, xaxs = "i", yaxs = "i", lwd = 2)
abline(v = mu)
axis(side = 1, at = Xlim)
axis(side = 1, mu, expression(mu))
text(Xlim[2], Ylim[2]*.8, adj = 1, paste("Population"))


### Samples
PlotSample(A[[first]], 1)
for(i in 2:nDatasets)
  PlotSample(A[[ii[i]]], i, col = myColors[i])

### Column 2: Sampling Distn, then 5 bootstrap samples
xd <- seq(from = XlimMean[1], to = XlimMean[2], length = 151)
plot(xd, dnorm(xd, sd = sd1/sqrt(n)),
     xlim = XlimMean, ylim = YlimMean, xlab = "", ylab = "",
     type = "l", axes = FALSE, main = "", yaxs = "i", lwd = 2)
axis(side = 1, XlimMean)
abline(v = mu, h = 0)
axis(side = 1, mu, expression(mu))
# Need axis call last, not abline, or arrows and text don't show up.
text(XlimMean[2], YlimMean[2]*.8, adj = 1, paste("Sampling\nDistribution"))
text(XlimMean[2], YlimMean[2]*.56, adj = 1, expression("of " * bar(x)))

for(i in 1:nDatasets){
  hist(bootstraps2[[i]], col = myColors[i],
       xlim = XlimMean, ylim = YlimMean, xlab = "", ylab = "",
       axes = FALSE, main = "", yaxs = "i", legend = FALSE)
  abline(v = mu)
  axis(side = 1, XlimMean)
  axis(side = 1, A1means[ii[i]], expression(bar(x)))
  if(bootstraps2[[i]]$observed > mu) {
    text(XlimMean[1], YlimMean[2]*.75, adj = 0,
         paste("Bootstrap\ndistribution\nfrom sample", i))
    text(XlimMean[1], YlimMean[2]*.44, adj = 0,
         expression("of " * bar(x)^"*"))
  } else {
    text(XlimMean[2], YlimMean[2]*.75, adj = 1,
         paste("Bootstrap\ndistribution\nfrom sample", i))
    text(XlimMean[2], YlimMean[2]*.44, adj = 1,
         expression("of " * bar(x)^"*"))
  }
}

### Column 3: 6 more bootstrap samples

for(i in 1:(1+nDatasets)){
  hist(bootstraps3[[i]], col = myColors[1],
       xlim = XlimMean, ylim = YlimMean, xlab = "", ylab = "",
       axes = FALSE, main = "", yaxs = "i", legend = FALSE)
  axis(side = 1, XlimMean)
  abline(v = mu)
  axis(side = 1, A1means[first], expression(bar(x)))
    text(XlimMean[2], YlimMean[2]*.7, adj = 1,
         paste0("Bootstrap\ndistribution\nfrom sample 1\n",
                "R = ", bootstraps3[[i]]$R))
}

#dev.off()
