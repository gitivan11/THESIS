# Produce figure with
#   column 1: population and samples
#   column 2: sampling distribution and bootstrap distributions
#   column 3: additional bootstrap distributions
# for the median, small odd n.

library(resample)
source("functions.R")

#--------------------------------------------------


rPopulation <- function(n){
  # Generate observations from population
  i <- sample(1:2, size = n, replace = TRUE)
  rnorm(n, mean = populationMean[i], sd = populationSd[i])
}
dPopulation <- function(x){
  # density from population
  (dnorm(x, populationMean[1], sd = populationSd[1]) +
   dnorm(x, populationMean[2], sd = populationSd[2]))/2
}
pPopulation <- function(x){
  # density from population
  (pnorm(x, populationMean[1], sd = populationSd[1]) +
   pnorm(x, populationMean[2], sd = populationSd[2]))/2
}
dMedian <- function(x, n){
  # Density for sampling distribution of the median
  # Number of observations from first half of population is Bi(n, .5)
  if(n %% 2 == 0) stop("n must be odd")
  n2 <- (n-1)/2
  a <- choose(n, n2) * (n2+1)
  a * pPopulation(x)^n2 * (1-pPopulation(x))^n2 * dPopulation(x)
}


myColors <- c(4, 2, 3, 5, 6) # blue first, then red green cyan magenta


jitter3 <-
function(x, spacing){
  # Perturb values of x if necessary, to maintain desired spacing.
  # Quick and dirty, perturbs more than necessary - lazy way to avoid an
  # infinite loop.
  #
  # Args:
  #   x       : sorted vector
  #   spacing : minimum spacing required in the output
  if(is.unsorted(x))
    x <- sort(x)

  # First modification (may get more complicated later)
  dx <- diff(x)
  if(all(dx >= spacing))
    return(x)
  change <- pmax(0, spacing-dx)
  x <- x - .5 * c(change, 0) + .5 * c(0, change)
  dx <- diff(x)
  if(all(dx >= spacing))
    return(x)

  # There were some instances where changing values created
  # new conflicts.  Resolve these by finding sets of observations
  # that are in conflict, and setting them to have the same mean
  # as originally, but with the desired spacing.
  N <- seq(along = x)
  repeat {
    j <- which(dx < spacing - 1e-6)[1]  # first obs that is in conflict
    M <- mean(x[j + 0:1]) - spacing/2   # where x[j] will be (maybe)
    xx <- (N - j) * spacing + M   # where obs would be with fixed spacing
    k <- which( (N <= j & x >= xx) | (N > j & x <= xx))
    x[k] <- xx[k] + (mean(x[k]) - mean(xx[k]))
    dx <- diff(x)
    if(all(dx >= spacing - 1e-6))
      return(x)
  }
}


#--------------------------------------------------
# Define population, limits, samples

#### Density for Population
# Population is 50% N(0, 1) and 50% N(4, 2^2); mu = 2, var = 4 + (1+4)/2 = 6.5
# Same as BootPix.R but shifted, and smaller n
# Similar to distributionGraph1, but less bimodal

n <- 15
populationMean <- c(0, 4)
populationSd <- c(1, 1.5)
mu <- mean(populationMean)
Xlim <- c(-3, 9)
median <- uniroot(function(x) pPopulation(x) - 0.5,
                  interval = c(Xlim[1], mu))$root


# Generate some random samples
set.seed(0)
A <- lapply(1:20, function(i, n) rPopulation(n), n = n)
A <- lapply(A, jitter3, spacing = -.05)
# Don't want any with values outside Xlim
bad <- (sapply(A, min) < Xlim[1]) | (sapply(A, max) > Xlim[2])
ii <- which(!bad)[1:5]
# 2:6

bootstraps2 <- lapply(ii,
                     function(i) bootstrap(A[[i]], median, seed = i, R = 10^4))
boot2summaries <- lapply(bootstraps2, function(x) table(x$replicates))
YlimBoot <- c(0, 1.04 * max(unlist(boot2summaries)))

smoothMedian <- function(x) {
  n <- length(x)
  median(x + rnorm(n, sd = sd(x)/sqrt(n)))
}
bootstraps3 <- lapply(ii,
                     function(i) bootstrap(A[[i]], smoothMedian, seed = i,
                                           R = 10^4))
for(i in 1:5) bootstraps3[[i]]$observed <- median(A[[ii[i]]])


PlotBars <- function(x, i, ...) {
  plot(x[1:2], c(0, 1.1), type = "n", axes = FALSE,
       xlab = "", xlim = Xlim, xaxs = "i", yaxs = "i", main = "")
  rect(xleft = x-.02, ybottom = 0, xright = x + .02, ytop = 1,
       col = myColors[i], ...,
       border = myColors[i], density = NA)
  text(Xlim[2], 1, adj = 1, paste("Sample", i))
  abline(v = median)
  axis(side = 1, at = Xlim)
  axis(side = 1, at = median(x), "m")
}

PlotBootBars <- function(i, ...) {
  xt <- boot2summaries[[i]]
  xv <- as.numeric(names(xt))
  plot(xv, xt, type = "n", axes = FALSE, ylim = YlimBoot,
       xlab = "", xlim = Xlim, xaxs = "i", main = "")
  axis(side = 1, at = Xlim)
  axis(side = 1, at = bootstraps2[[i]]$observed, "m")
  rect(xleft = xv-.02, ybottom = 0, xright = xv + .02, ytop = xt,
       col = myColors[i], ..., border = myColors[i], density = NA)
  text(Xlim[2], YlimBoot[2]*.8, adj = 1,
       paste("Bootstrap\ndistribution\nfrom sample", i))
  abline(v = median)
}


#--------------------------------------------------
# Make distributionGraph3

# Either 5 or 4 datasets (below the population)
nDatasets <- 4  # or 5


#pdf("../figures/distributionGraph3.pdf", height = 1+nDatasets, width = 6)
par(mfcol = c(1+nDatasets, 3), cex = .5, mar = c(4.1, .5, .1, .5))

### Column 1: population density, then 5 samples
x <- seq(from = Xlim[1], to = Xlim[2], length = 151)
plot(x, dPopulation(x), xlab = "", ylab = "",
     type = "l", xlim = Xlim, ylim = range(dPopulation(x))*1.04,
     axes = FALSE, xaxs = "i", yaxs = "i", lwd = 2)
abline(v = median)
axis(side = 1, at = Xlim)
axis(side = 1, median, "M") # population median
text(Xlim[2], par("usr")[4] * .8, adj = 1, "Population\nmedian = M")


### Samples
for(i in 1:nDatasets)
  PlotBars(A[[ii[i]]], i)


### Column 2: Sampling Distn, then 5 bootstrap samples
y <- dMedian(x, n)
Ylim2 <- c(0, max(y)*1.04)
plot(x, y,
     xlim = Xlim, ylim = Ylim2, xlab = "", ylab = "",
     type = "l", axes = FALSE, main = "", yaxs = "i", lwd = 2)
abline(v = median)
axis(side = 1, at = Xlim)
axis(side = 1, median, "M") # population median
text(Xlim[2], par("usr")[4] * .7, adj = 1,
     paste("Sampling\ndistribution\nof median m\nwith n =", n))

for(i in 1:nDatasets)
  PlotBootBars(i)

### Column 3: Smoothed bootstrap distributions
plot(0:1, 0:1, type = "n", axes = FALSE, xlab = "", ylab = "", main = "")
text(.5, .5, adj = .5, "Smoothed Bootstrap")

for(i in 1:nDatasets){
  hist(bootstraps3[[i]], col = myColors[i],
       xlim = Xlim, ylim = Ylim2 * 1.3, xlab = "", ylab = "",
       axes = FALSE, main = "", yaxs = "i", legend = FALSE)
  abline(v = median)
  axis(side = 1, Xlim)
  axis(side = 1, bootstraps2[[i]]$observed, "m")
  text(Xlim[2], Ylim2[2]*.75, adj = 1,
       paste("Smoothed\nbootstrap\ndistribution\nfrom sample", i))
}

#dev.off()
