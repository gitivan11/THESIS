# Bootstrap tilting.
# Using Verizon CLEC example. Also Verizon ILEC for comparison.

library(resample)
library(aggregate)  # From http://www.timhesterberg.net/r-packages
source("functions.R")

data(Verizon) # This comes with library(resample)

CLEC <- sort(with(Verizon, Time[Group == "CLEC"]))
ILEC <- sort(with(Verizon, Time[Group == "ILEC"]))
n.CLEC <- length(CLEC)
n.ILEC <- length(ILEC)

# w = weights on observations, 1 to n. Sum to 1.
# W = weights on bootstrap samples.


weight <- function(tau, x) {
  # Tilting weights
  w <- 1 / (1 - tau * (x - mean(x)))
  stopifnot(all(w > 0))
  w / sum(w)
}

# I plan to show 5 distributions.
# The two most extreme should be endpoints of a 95% tilting CI.
# The middle one should be no tilting; that would be roughly a one-sided 50% CI,
# though not exactly. For the 2nd and 4th, interpolate.

#bootstrap(CLEC, mean, seed = 0)
#bootstrap(ILEC, mean, seed = 1)
# Gives same results. But don't define, because elsewhere CLEC isn't sorted.

# Do bootstrapping by hand, in order to save replicates
r <- 10^4
set.seed(0)
CLEC.indices <- samp.bootstrap(23, r)
set.seed(1)
ILEC.indices <- samp.bootstrap(1664, r)
CLEC.means <- indexMeans(CLEC, CLEC.indices)
ILEC.means <- indexMeans(ILEC, ILEC.indices)

RelativeLikelihood <- function(indices, w) {
  # Calculate the relative likelihood of bootstrap samples under w and w0
  n <- nrow(indices)
  # Importance sampling weight for each bootstrap sample
  exp(indexSums(log(n * w), indices))
}
plot(CLEC.means, RelativeLikelihood(CLEC.indices, weight(0, CLEC)))
plot(CLEC.means, RelativeLikelihood(CLEC.indices, weight(-.001, CLEC)))
# 1.3 to .6
plot(CLEC.means, RelativeLikelihood(CLEC.indices, weight(.001, CLEC)))
# .8 to 1.6

# I'll use the "integration" importance sampling estimate here
# Find the (positive) tau with P(Xbar^* < xbar) = .025
# Find the (negative) tau with P(Xbar^* > xbar) = .025

TailProb <- function(tau, x, indices, lower.tail = TRUE, subtract = 0){
  # Probability of Xbar^* < xbar, when sampling using tau
  observed <- mean(x)
  statistics <- indexMeans(x, indices)
  use <- ifelse1(lower.tail, statistics < observed, statistics > observed)
  prob <- RelativeLikelihood(indices, weight(tau, x))
  sum(prob[use]) / ncol(indices) - subtract
}
1 / range(CLEC - mean(CLEC)) # -0.06057254  < tau < 0.01252962
1 / range(ILEC - mean(ILEC)) # -0.118883297 < tau < 0.005458861
TailProb(0, CLEC, CLEC.indices) # .5438
TailProb(0, CLEC, CLEC.indices, lower.tail = FALSE) # .4561
TailProb(0.001, CLEC, CLEC.indices) # .505
TailProb(0.002, CLEC, CLEC.indices) # .463
TailProb(0.003, CLEC, CLEC.indices) # .417
TailProb(0.004, CLEC, CLEC.indices) # .367
TailProb(0.005, CLEC, CLEC.indices) # .312
TailProb(0.006, CLEC, CLEC.indices) # .253
TailProb(0.007, CLEC, CLEC.indices) # .191
TailProb(0.008, CLEC, CLEC.indices) # .127
TailProb(0.009, CLEC, CLEC.indices) # .068
TailProb(0.010, CLEC, CLEC.indices) # .023
# The tau to get to .025 is surprisingly close to the constraint
TailProb(-0.01, CLEC, CLEC.indices, lower.tail = FALSE) # .20
TailProb(-0.02, CLEC, CLEC.indices, lower.tail = FALSE) # .089
TailProb(-0.03, CLEC, CLEC.indices, lower.tail = FALSE) # .032
TailProb(-0.04, CLEC, CLEC.indices, lower.tail = FALSE) # .008

tau.CLEC1 <- uniroot(TailProb, c(-.04, -.03), x = CLEC, indices = CLEC.indices,
                     subtract = .025, lower.tail = FALSE)$root
tau.CLEC2 <- uniroot(TailProb, c(.009, .010), x = CLEC, indices = CLEC.indices,
                     subtract = .025)$root
tau.CLEC <- c(tau.CLEC1 * c(1, .5), 0, tau.CLEC2 * c(.5, 1))
tau.CLEC
# -0.032110920 -0.016055460  0.000000000  0.004969482  0.009938965
CLEC.p <- c(.975, NA, mean(CLEC.means < mean(CLEC)), NA, .025)
CLEC.p[2] <- TailProb(tau.CLEC[2], CLEC, CLEC.indices)
CLEC.p[4] <- 1 - TailProb(tau.CLEC[4], CLEC, CLEC.indices, lower.tail = FALSE)
plot(tau.CLEC, CLEC.p)

tau.ILEC1 <- uniroot(TailProb, c(-.1, 0), x = ILEC, indices = ILEC.indices,
                     subtract = .025, lower.tail = FALSE)$root
tau.ILEC2 <- uniroot(TailProb, c(0, .0054), x = ILEC, indices = ILEC.indices,
                     subtract = .025)$root
ILEC.p <- c(.975, NA, mean(ILEC.means < mean(ILEC)), NA, .025)
ILEC.p[2] <- TailProb(tau.ILEC[2], ILEC, ILEC.indices)
ILEC.p[4] <- 1 - TailProb(tau.ILEC[4], ILEC, ILEC.indices, lower.tail = FALSE)

tau.ILEC <- c(tau.ILEC1 * c(1, .5), 0, tau.ILEC2 * c(.5, 1))
tau.ILEC
# -0.003656289 -0.001828144  0.000000000  0.001403094  0.002806189
plot(tau.ILEC, ILEC.p)


AdjustProbabilities <- function(W, method = "regression") {
  # W is a vector of probabilities; adjust to mean 1
  Wbar <- mean(W)
  if(Wbar == 1)
    return(W)
  n <- length(W)
  if(method == "ratio")
    return(W / Wbar)
  Wd <- W - Wbar
  if(method == "regression") {
    return(W * (1 + n * (1-Wbar) * Wd / sum(Wd^2)))
  }
  if(method == "ML") {
    # Maximum likelihood. maximize p s.t. sum(p*W) = 1, sum(p)=1, p>=0
    # p = c / (1 - tau * (W - Wbar))
    # Then return p*W*n
    f <- function(tau, W, Wd) { # function to solve for tau
      p <- 1 / (1 - tau * Wd) # unnormalized
      sum(p * W) / sum(p) - 1
    }
    tau <- uniroot(f, range((2*n+1)/(2*n+2) / range(Wd)), W = W, Wd = Wd)$root
    p <- 1 / (1 - tau * Wd)
    return(p * W / mean(p))
  }

  # Another method:
  # Primarily adjust the largest values.
  return(NA)
}

CLEC.w <- matrix(NA, n.CLEC, 5)
ILEC.w <- matrix(NA, n.ILEC, 5)
CLEC.W <- matrix(NA, ncol(CLEC.indices), 5)
ILEC.W <- matrix(NA, ncol(ILEC.indices), 5)
for(i in 1:5) {
  CLEC.w[, i] <- weight(tau.CLEC[i], CLEC)
  ILEC.w[, i] <- weight(tau.ILEC[i], ILEC)
  CLEC.W[, i] <- RelativeLikelihood(CLEC.indices, CLEC.w[, i])
  ILEC.W[, i] <- RelativeLikelihood(ILEC.indices, ILEC.w[, i])
}
colSums(CLEC.w)
colSums(ILEC.w)
# all 1's (good)
CLEC.tiltMeans <- colSums(CLEC * CLEC.w)
ILEC.tiltMeans <- colSums(ILEC * ILEC.w)
CLEC.tiltMeans
# 10.65714 13.04979 16.50913 19.16571 28.87603
ILEC.tiltMeans
# 7.761783 8.057521 8.411611 8.748604 9.196132

colMeans(CLEC.W) # 1.0333450 1.0032266 1.0000000 0.9996052 0.9031924
colMeans(ILEC.W) # 0.9924466 0.9980302 1.0000000 1.0073979 0.9553067
# Need to adjust. Before I sorted CLEC, one of the means was 1.4, with max=6221
CLEC.Wreg <- matrix(NA, ncol(CLEC.indices), 5)
ILEC.Wreg <- matrix(NA, ncol(ILEC.indices), 5)
for(i in 1:5) {
  CLEC.Wreg[, i] <- AdjustProbabilities(CLEC.W[, i], "regression")
  ILEC.Wreg[, i] <- AdjustProbabilities(ILEC.W[, i], "regression")
}
colMeans(CLEC.Wreg)
colMeans(ILEC.Wreg)
devAskNewPage(ask = TRUE)
for(i in 1:5) {
  PlotSort(CLEC.W[, i], main = paste("CLEC column", i, "mean", round(mean(CLEC.W[, i]), 3)))
  PointsSort(CLEC.Wreg[, i], col = 2, pch=".")
  plot(sqrt(CLEC.W[, i]), sqrt(CLEC.Wreg[, i]), pch="."); abline(0, 1)
  PlotSort(ILEC.W[, i], main = paste("ILEC column", i, "mean", round(mean(ILEC.W[, i]), 3)))
  PointsSort(ILEC.Wreg[, i], col = 2, pch=".")
  plot(sqrt(ILEC.W[, i]), sqrt(ILEC.Wreg[, i]), pch="."); abline(0, 1)
}

CLEC.Wml <- matrix(NA, ncol(CLEC.indices), 5)
ILEC.Wml <- matrix(NA, ncol(ILEC.indices), 5)
for(i in 1:5) {
  CLEC.Wml[, i] <- AdjustProbabilities(CLEC.W[, i], "ML")
  ILEC.Wml[, i] <- AdjustProbabilities(ILEC.W[, i], "ML")
}
colMeans(CLEC.Wml)
colMeans(ILEC.Wml)
# Not exactly 1, but close
for(i in 1:5) {
  PlotSort(CLEC.W[, i], main = paste("CLEC column", i, "mean", round(mean(CLEC.W[, i]), 3)))
  PointsSort(CLEC.Wreg[, i], col = 2, pch=".")
  PointsSort(CLEC.Wml[, i], col = 3, pch=".")
  plot(sqrt(CLEC.W[, i]), sqrt(CLEC.Wreg[, i]), col=2, pch=".")
  points(sqrt(CLEC.W[, i]), sqrt(CLEC.Wml[, i]), col=3, pch=".")
  abline(0, 1)
  PlotSort(ILEC.W[, i], main = paste("ILEC column", i, "mean", round(mean(ILEC.W[, i]), 3)))
  PointsSort(ILEC.Wreg[, i], col = 2, pch=".")
  PointsSort(ILEC.Wml[, i], col = 3, pch=".")
  plot(sqrt(ILEC.W[, i]), sqrt(ILEC.Wreg[, i]), col=2, pch=".")
  points(sqrt(ILEC.W[, i]), sqrt(ILEC.Wml[, i]), col=3, pch=".")
  abline(0, 1)
}
# Regression and ML adjustments are similar. No big advantage for ML.
# Use Regression, seems safer.
devAskNewPage(ask = FALSE)

# Actually do weighted sampling
set.seed(1)
CLEC.means2 <- matrix(NA, r, 5)
ILEC.means2 <- matrix(NA, r, 5)
for(i in 1:5) {
  temp <- matrix(sample(n.CLEC, size = r * n.CLEC,
                        replace = TRUE, prob = CLEC.w[,i]),
                 n.CLEC, r)
  CLEC.means2[, i] <- indexMeans(CLEC, temp)
  temp <- matrix(sample(n.ILEC, size = r * n.ILEC,
                        replace = TRUE, prob = ILEC.w[,i]),
                 n.ILEC, r)
  ILEC.means2[, i] <- indexMeans(ILEC, temp)
}

# Double-check the one-sided P-values; the two extreme should be 2.5%
colMeans(CLEC.means2 < mean(CLEC))
# 0.9764 0.8721 0.5345 0.3119 0.0246
colMeans(ILEC.means2 < mean(ILEC))
# 0.9732 0.8516 0.5054 0.1820 0.0251


# p.CLEC <- mean(CLEC.means < mean(CLEC)) # .5438
# p.ILEC <- mean(ILEC.means < mean(ILEC)) # .5087

# weight(tau, x) # returns n-vector of weights
# RelativeLikelihood(indices, w) # return r-vector of relative likelihood
# TailProb(tau, x, indices, lower.tail = TRUE, subtract = 0) # boot CDF
# AdjustProbabilities(W, method = "regression")

#--------------------------------------------------
# Plots

cumsum2 <- function(x) cumsum(x) - x/2
o.CLEC <- order(CLEC.means)
o.ILEC <- order(ILEC.means)
colorVec <- c("blue", "cyan", "gray", "magenta", "red")
colorVec <- c("blue", "lightblue", "gray", "pink", "red")
colorVec <- c("blue", "lightblue", "gray", "magenta", "red")

# # Bootstrap distributions using impsamp. Use Wreg instead of W (W blows up)
# o <- o.CLEC
#
# #par(mfcol = c(6, 3), mex = .5, cex = .5)
# for(i in 1:5) {
#   plot(density(CLEC.means, weights = CLEC.Wreg[, i] / r,
#                bw = sd(CLEC.means)/6), type = "l",
#        main = "", xlab = "") #"CLEC bootstrap means")
#   AddLines(i)
# }
# # At bottom, 5 CDFs superimposed, xbar* vs cumsum(W)
# plot(range(CLEC.means), 0:1, type = "n", xlab = "CLEC bootstrap means",
#      ylab = "CDF")
# for(i in 1:5)
#   lines(CLEC.means[o], cumsum2(CLEC.Wreg[o, i]) / r,
#         col = colorVec[i])
# abline(v = mean(CLEC), col = "gray")




#--------------------------------------------------
# pdf("../figures/TiltCLEC.pdf", height = 8, width = 6)
par(mfcol = c(6, 3), mex = .5, cex = .5)
par(oma = c(1, 1, 0, 0), mar = c(4.6, 4.6, 1.1, 1.1))
# 3 columns: x vs w, density of data, density of bootstrap means
# At bottom of each column show a plot that summarizes

AddLines <- function(i) {
  axis(side = 1, at = sum(CLEC.w[, i] * CLEC), expression(mu),
       col = colorVec[i])
  abline(v = mean(CLEC), col = "gray")
  abline(v = CLEC.tiltMeans[i], col = colorVec[i])
}

# 5 individual scatterplots, x vs w
for(i in 1:5) {
  plot(CLEC, CLEC.w[, i], ylim = range(0, CLEC.w), type = "n",
       ylab = "Probability", xlab = "")
  abline(h = 1/n.CLEC, col = "gray")
  points(CLEC, CLEC.w[, i])
  segments(CLEC, CLEC.w[, i], y1 = 0, ylim = range(CLEC.w))
  AddLines(i)
}
# at bottom, matplot
matplot(CLEC, CLEC.w, type = "b", lty=1, pch=1,
        col = colorVec,
        ylim = range(0, CLEC.w), xlab = "CLEC (hours)", ylab = "Probability")


#par(mfcol = c(6, 3), mex = .5, cex = .5)
# 5 density plots, weighted x
temp <- list()
for(i in 1:5)
  temp[[i]] <- density(CLEC, weights = CLEC.w[, i], bw = sd(CLEC)/sqrt(n.CLEC))
temp2 <- max(sapply(temp, function(x) max(x$y)))
for(i in 1:5) {
  plot(temp[[i]], main="", xlab = "", ylim = c(0, temp2))
  AddLines(i)
}
# At bottom, 5 CDFs superimposed, x vs cumsum(w)
plot(range(temp[[1]]$x), 0:1, type = "n", xlab = "CLEC (hours)", ylab = "CDF")
for(i in 1:5)
  lines(CLEC, cumsum2(CLEC.w[, i]), col = colorVec[i])


# Bootstrap distributions using weighted sampling
temp <- list()
for(i in 1:5)
  temp[[i]] <- density(CLEC.means2[, i], bw = sd(CLEC.means)/6)
temp2 <- max(sapply(temp, function(x) max(x$y)))
for(i in 1:5) {
  plot(temp[[i]], ylim = c(0, temp2),
       type = "l", main = "", xlab = "", xlim = range(CLEC.means))
  AddLines(i)
  if(i == 1)
    legend("topright", lty=rep(1,2), col = colorVec[c(1,3)],
       legend = expression(mu, bar(x)))
}
# At bottom, 5 CDFs superimposed, xbar* vs cumsum(W)
plot(range(CLEC.means), 0:1, type = "n", xlab = "CLEC bootstrap means",
     ylab = "CDF")
for(i in 1:5)
  lines(sort(CLEC.means2[, i]), ppoints(r), col = colorVec[i])
abline(v = mean(CLEC), col = "gray")

# dev.off()

#--------------------------------------------------
# pdf("../figures/TiltILEC.pdf", height = 8, width = 6)
par(mfcol = c(6, 3), mex = .5, cex = .5)
par(oma = c(1, 1, 0, 0), mar = c(4.6, 4.6, 1.1, 1.1))
# 3 columns: x vs w, density of data, density of bootstrap means
# At bottom of each column show a plot that summarizes

AddLines <- function(i) {
  axis(side = 1, at = sum(ILEC.w[, i] * ILEC), expression(mu),
       col = colorVec[i])
  abline(v = mean(ILEC), col = "gray")
  abline(v = ILEC.tiltMeans[i], col = colorVec[i])
}

# 5 individual scatterplots, x vs w
for(i in 1:5) {
  plot(ILEC, ILEC.w[, i], ylim = range(0, ILEC.w), type = "n",
       ylab = "Probability", xlab = "")
  abline(h = 1/n.ILEC, col = "gray")
  points(ILEC, ILEC.w[, i])
  segments(ILEC, ILEC.w[, i], y1 = 0, ylim = range(ILEC.w))
  AddLines(i)
}
# at bottom, matplot
matplot(ILEC, ILEC.w, type = "b", lty=1, pch=1,
        col = colorVec,
        ylim = range(0, ILEC.w), xlab = "ILEC (hours)", ylab = "Probability")


#par(mfcol = c(6, 3), mex = .5, cex = .5)
# 5 density plots, weighted x
temp <- list()
for(i in 1:5)
  temp[[i]] <- density(ILEC, weights = ILEC.w[, i], bw = sd(ILEC)/sqrt(n.ILEC))
temp2 <- max(sapply(temp, function(x) max(x$y)))
for(i in 1:5) {
  plot(temp[[i]], main="", xlab = "", ylim = c(0, temp2))
  AddLines(i)
}
# At bottom, 5 CDFs superimposed, x vs cumsum(w)
plot(range(temp[[1]]$x), 0:1, type = "n", xlab = "ILEC (hours)", ylab = "CDF")
for(i in 1:5)
  lines(ILEC, cumsum2(ILEC.w[, i]), col = colorVec[i])


# Bootstrap distributions using weighted sampling
temp <- list()
for(i in 1:5)
  temp[[i]] <- density(ILEC.means2[, i], bw = sd(ILEC.means)/6)
temp2 <- max(sapply(temp, function(x) max(x$y)))
for(i in 1:5) {
  plot(temp[[i]], ylim = c(0, temp2),
       type = "l", main = "", xlab = "", xlim = range(ILEC.means))
  AddLines(i)
  if(i == 1)
    legend("topright", lty=rep(1,2), col = colorVec[c(1,3)],
       legend = expression(mu, bar(x)))
}
# At bottom, 5 CDFs superimposed, xbar* vs cumsum(W)
plot(range(ILEC.means), 0:1, type = "n", xlab = "ILEC bootstrap means",
     ylab = "CDF")
for(i in 1:5)
  lines(sort(ILEC.means2[, i]), ppoints(r), col = colorVec[i])
abline(v = mean(ILEC), col = "gray")

# dev.off()

#--------------------------------------------------
# Mean-variance relationship
