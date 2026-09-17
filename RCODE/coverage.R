# Coverage of confidence intervals
# This file was used for the arXiv version, and the first TAS version

library(Hmisc)
library(aggregate)  # From http://www.timhesterberg.net/r-packages

# TODO: add BCa and ABC intervals

#--------------------------------------------------
nn <- c(5, 10, 20, 40, 80)
nn1 <- c(5, 10, 20, 40, 80, 160)
nn2 <- c(100, 200, 400, 1000)
nn3 <- c(5, 10, 20, 40, nn2)
nn4 <- c(2000, 4000)
nn5 <- c(nn3, nn4)
nn6 <- c(6000, 8000)
nn7 <- c(10000, 15000, 20000)


#--------------------------------------------------
# Narrowness factors

# I expect these one-sided miss probabilities for z interval:
pt(qnorm(.025), nn - 1)
# 0.06077983 0.04082457 0.03241674

# For percentile interval, ignoring random skewness, expect
pt(qnorm(.025) * sqrt((nn-1)/nn), nn - 1)
# 0.07723296 0.04795364 0.03565071

# My adjusted percentile interval uses:
probs2 <- pnorm(qt(.025, nn-1) * sqrt(nn/(nn-1)))
probs2
# [1] 0.0009541006 0.0085506390 0.0158808299 0.0202575364 0.0225881945
pt(qnorm(probs2) * sqrt((nn-1)/nn), nn - 1)
# All .025 (actual coverage will be worse, due to skewness variation,
# and extra variation due to estimating width using quantiles)

# table:narrowness
latex(round(cbind(nn,
                  sqrt((nn-1)/nn),
                  qnorm(.975) / qt(.975, nn-1)),
            3), file = "")
# Add another column, with an extra digit
latex(cbind(round(cbind(nn,
                        sqrt((nn-1)/nn),
                        qnorm(.975) / qt(.975, nn-1)), 3),
            round(pnorm(sqrt(nn/(nn-1)) * qt(.025, nn-1)), 4)),
      file = "")
# Insert column with actual non-coverage xbar +- z sigmahat/sqrt(n)
latex(cbind(round(cbind(nn,
                        sqrt((nn-1)/nn),
                        qnorm(.975) / qt(.975, nn-1),
                        pt(qnorm(.025)*sqrt((nn-1)/nn), nn-1)
                        ),
                  3),
            round(pnorm(sqrt(nn/(nn-1)) * qt(.025, nn-1)), 4)),
      file = "")



#--------------------------------------------------

meanSE <- function(x) {
  # Mean and standard error for the mean
  if(is.matrix(x))
    return(rbind(colMeans(x), colStdevs(x)/sqrt(nrow(x))))
  c(mean(x), sd(x)/sqrt(length(x)))
}

ComputeIntervals <- function(x, xbar = mean(x), n = length(x)) {
  # Compute a number of intervals
  # x is a single sample.
  # xbar can be supplied if known.
  alpha <- c(.025, .975)
  talpha <- qt(alpha, n-1)
  probs2 <- pnorm(talpha * sqrt(n/(n-1))) # adjusted percentile
  s <- sd(x)
  sqrtn <- sqrt(n)
  ciT <- xbar + talpha * s/sqrtn
  A <- mean((x-xbar)^3) / s^3 / 6 / sqrtn
  ciJohnson <- xbar + s/sqrtn * (talpha + A * (2*talpha^2+1))
  # Bootstrap percentile endpoints
  m <- matrix(sample(x, size = R*n, replace = TRUE), n)
  bootMeans <- colMeans(m)
  bootS <- colStdevs(m)
  ciTBoot <- xbar + talpha * sd(bootMeans)
  ciPercentile <- quantile(bootMeans, probs = alpha, type = 6)
  ciPercentileRev <- 2*xbar - rev(ciPercentile)
  ciPercentileAdj <- quantile(bootMeans, probs = probs2, type = 6)
  ciBoott <- xbar - quantile((bootMeans-xbar)/bootS, probs = rev(alpha),
                             type = 6) * s
  rbind(t = ciT,
        tSkew = ciJohnson,
        tBoot = ciTBoot,
        perc = ciPercentile,
        reverse = ciPercentileRev,
        expanded = ciPercentileAdj, # This is called percAdj in some saved objects
        bootT = ciBoott)
}

# TODO: redo the figures to rename percAdj to percExp or expanded

PlotResults <- function(results, ylim = NULL, ..., average = FALSE,
                        legend = TRUE,
                        ylab = "One-sided Non-coverage", xlab = "n",
                        legend.cex = par("cex"),
                        scale = c("inverseRoot", "inverse", "inverse4Root"),
                        nMax = max(n), nLim = range(n), axis2 = TRUE,
                        omit = NULL, ltyVec = 1:k) {
  # Plot confidence interval coverage (both one-sided coverages)
  # Args:
  #   results: e.g. gammaResults or normalResults
  #   average: if TRUE average the one-sided non-coverages and only plot the avg
  #   scale: The x axis is plotted on a transformed scale:
  #     "inverseRoot" -> -1/sqrt(n)
  #     "inverse" -> -1/n
  #     "inverse4Root" -> -1/n^(1/4)
  #   nMax: For n larger than this, only lines are used, not points
  #   omit: e.g. 'tSkew' to omit that one.
  #   ltyVec: lty for: t S B p r e T  (by default, t and T are both 1)
  #
  # Details:
  #   Plotting is done on a transformed scale
  scale <- match.arg(scale)
  Transform <- switch(scale,
                      inverseRoot = function(n) -1/sqrt(n),
                      inverse = function(n) -1/n,
                      inverse4Root = function(n) -1/n^.25)

  means <- lapply(results, colMeans)
  means1 <- do.call("rbind", lapply(means, function(x) x[, 1]))
  means2 <- 1 - do.call("rbind", lapply(means, function(x) x[, 2]))
  n <- as.numeric(names(results))
  o <- order(n)
  n <- n[o]
  means1 <- means1[o, ]
  means2 <- means2[o, ]
  if(is.null(ylim))
    ylim <- range(means1, means2)
  k <- ncol(means1)
  x <- Transform(n)
  if(average)
    means2 <- (means2+means1)/2
  #colors <- c(1:6, 8) # avoid yellow, color 7
  colors <- rep(1, 7)
  pchVec <- c("t", "S", "B", "p", "r", "e", "T")
  showPoints <- (n <= nMax)
  if(is.na(nLim[2]))
    nLim[2] <- max(n)
  keep <- !(colnames(means2) %in% omit)
  matplot(x[showPoints], means2[showPoints, keep , drop = FALSE],
          pch = pchVec[keep],
          lty = ltyVec, col = colors[keep],
          xlim = range(Transform(nLim)), ylim = ylim, type = "b", axes = FALSE,
          xlab = xlab, ylab = ylab, ...)
  if(!all(showPoints)) {
    showPoints[n == nMax] <- FALSE
    matlines(x[!showPoints], means2[!showPoints, keep, drop = FALSE],
             lty = ltyVec, col = colors[keep])
  }
  if(axis2)
    axis(side = 2)
  axis(side = 1, at = Transform(n), labels = as.character(n))
  if(!average)
    matlines(x, means1[, keep, drop = FALSE], lty = ltyVec,
                       col = colors[keep], type = "l")
  if(axis2)
    axis(side = 4, at = 0.025)
  if(legend)
    legend("topright", lty = ltyVec, col = colors[keep],
           legend = paste0(pchVec[keep], ": ", colnames(means1)[keep]),
           cex = legend.cex)
  # legend("top", legend = "Upper endpoints have numbers")
  invisible(NULL)
}

PlotResults(normalResults, ylim = c(0, 0.10), nMax = 40)
PlotResults(normalResults, ylim = c(0, 0.10), scale = "inverse", nMax = 40)
PlotResults(save.gammaResults6, ylim = c(0, 0.10))
PlotResults(save.gammaResults6, ylim = c(0, 0.10), nMax = 200)
PlotResults(save.gammaResults6, ylim = c(0, 0.10), scale = "inverse")

#--------------------------------------------------
# Confidence interval coverage for small n, normal population

# Variance Reduction Trick: xbar and (x - xbar) are independent.
# Draw samples of v = (x-xbar), then integrate over Xbar
# A confidence interval for (xbar + v) = xbar + (interval based on v)

# Given any confidence interval based on v, compute expected value wrt Xbar:
#   P(endpoint + Xbar > theta) # where endpoint is based on v
#     = P(Xbar > theta - endpoint)

# In the inner part of the loop compute the endpoints.
# Outside, compute pnorm

# CDF of Xbar
FF <- function(x, n) pnorm(x, sd = sqrt(1/n))
mu <- 0

replications <- 10^4 # 10^2 for testing, 10^4 for real
R <- 10000 # bootstrap replications
CItemplate <- ComputeIntervals((1:10)^2)


normalCI <- list()
normalResults <- list()
normalIResults <- list()

# Create a 3-d array (replications, intervals, 2 endpoints
CIarray <- array(dim = c(replications, dim(CItemplate)),
                 dimnames = c(list(NULL), dimnames(CItemplate)))

for(n in nn) {
  cat("\nn =", n, "  ")
  set.seed(0)
  save.xbar <- rep(NA, replications)
  for(i in 1:replications) {
    if(i %in% (1:10 * (replications/10))) cat(i, " ")
    x <- rnorm(n)
    v <- x - mean(x) # Compute intervals based on v, not x.
    CIarray[i, , ] <- ComputeIntervals(v, xbar = 0, n)
    save.xbar[i] <- mean(x)
  }
  # Compute E(endpoint < mu)
  normalCI[[as.character(n)]] <- CIarray
  normalIResults[[as.character(n)]] <- (CIarray + save.xbar > mu)
  normalResults[[as.character(n)]] <- 1-FF(-CIarray, n) # mu - CIarray
}
colMeans(normalIResults[[1]])
colMeans(normalResults[[1]])
lapply(normalIResults, colMeans)
lapply(normalResults, colMeans)

# Use this for standard errors
lapply(normalResults, colStdevs)

# For normal results, average the non-coverage on the two sides
normalResults1 <- lapply(normalResults, function(x) (x[, ,1]+1-x[, ,2])/2)
lapply(normalResults1, colMeans)

PlotResults(normalIResults, ylim = c(0, 0.10))
PlotResults(normalResults, average = TRUE, ylim = c(.02, .05))



if(FALSE){
  save.normalResults6 <- normalResults
  save.normalIResults6 <- normalIResults
  save(list = c("replications", "R", "save.normalResults6"),
       file = "normalResults6.RData")
}

# save.*4 has R=10^4, 10^4 replications, and nn1
# save.*5 has R=10^4, 10^4 replications, and nn2
# save.*6 has R=10^4, 10^4 replications, nn5, and wider quantiles



# latex(round(... , 4), file = "")

# > lapply(normalResults, colMeans)
# $`5`
#               2.5%     97.5%
# t       0.02586050 0.9741395
# tAdj    0.03113588 0.9695317
# perc    0.08268674 0.9175258
# reverse 0.08247423 0.9173133
# percAdj 0.04083073 0.9598213
# boott   0.02680735 0.9742070
#
# $`10`
#               2.5%     97.5%
# t       0.02533283 0.9746672
# tAdj    0.02821614 0.9716871
# perc    0.04957614 0.9503487
# reverse 0.04965128 0.9504239
# percAdj 0.02733136 0.9726178
# boott   0.02551210 0.9745296
#
# $`20`
#               2.5%     97.5%
# t       0.02519000 0.9748100
# tAdj    0.02633731 0.9737130
# perc    0.03627467 0.9637253
# reverse 0.03627468 0.9637253
# percAdj 0.02570189 0.9742847
# boott   0.02530394 0.9747589
#
# $`40`
#               2.5%     97.5%
# t       0.02483949 0.9751605
# tAdj    0.02518403 0.9747887
# perc    0.03012186 0.9698821
# reverse 0.03011786 0.9698781
# percAdj 0.02504235 0.9749632
# boott   0.02490153 0.9750556




# Older

# 10^4 replications with R=1000
round(data.frame(lapply(nn, f, results = save.normalResults1)), 4)
#          X5    X10    X20    X40    X80
# t    0.0249 0.0246 0.0249 0.0248 0.0250
# z    0.0608 0.0403 0.0323 0.0283 0.0268
# tb   0.0340 0.0299 0.0277 0.0263 0.0258
# p.1  0.0822 0.0495 0.0368 0.0309 0.0285
# p.2  0.0821 0.0494 0.0368 0.0309 0.0284
# a.1  0.0419 0.0277 0.0263 0.0259 0.0260
# a.2  0.0419 0.0277 0.0262 0.0258 0.0259
# bt.1 0.0263 0.0257 0.0257 0.0257 0.0259
# bt.2 0.0261 0.0257 0.0258 0.0256 0.0259

# 10^4 replications with R=10000
round(data.frame(lapply(nn, f, results = save.normalResults2)), 4)
#          X5    X10    X20    X40    X80
# t    0.0259 0.0253 0.0252 0.0248 0.0250
# z    0.0622 0.0413 0.0326 0.0284 0.0267
# tb   0.0350 0.0306 0.0280 0.0263 0.0257
# p.1  0.0827 0.0496 0.0363 0.0301 0.0276
# p.2  0.0825 0.0497 0.0363 0.0301 0.0276
# a.1  0.0408 0.0273 0.0257 0.0250 0.0251
# a.2  0.0402 0.0274 0.0257 0.0250 0.0251
# bt.1 0.0268 0.0255 0.0253 0.0249 0.0251
# bt.2 0.0258 0.0255 0.0252 0.0249 0.0250

# Estimate the effective number of replications, with the variance reduction

temp.p <- do.call("rbind", lapply(save.normalResults4, colMeans))
temp.var <- do.call("rbind", lapply(save.normalResults4, colVars))

# Multiplier
temp.p * (1-temp.p) / temp.var
# 9.6 for n=5, over 500 for n=160

# Effective number of replications
temp.p * (1-temp.p) / (temp.var / replications)
# at least 96K replications

#--------------------------------------------------
# Exponential population

# Variance Reduction Trick: xbar and (x / xbar) are independent.
# Draw samples of v = (x/xbar), then integrate over xbar
# A confidence interval for (xbar * v) = xbar * (interval based on v)

# Given any confidence interval based on v, compute expected value wrt Xbar:
#   P(endpoint * Xbar > theta) # where endpoint is based on v
#     = P(Xbar > theta / endpoint)

# In the inner part of the loop compute the endpoints.
# Outside, compute pgamma


# CDF of Xbar
FG <- function(x, n) pgamma(x, shape = n, rate = n)
mu <- 1

replications <- 10^4 # 10^2 for testing, 10^4 for real, 10^3 for n=6K & 8K
R <- 10000 # bootstrap replications

if(FALSE){
  nn <- nn6
  replications <- 10^3
}
if(FALSE){
  nn <- nn7
  replications <- .5*10^3
}
if(FALSE){
  nn <- nn6
  replications <- 5e3
}


gammaCI <- list()
gammaResults <- list()
gammaIResults <- list()

# Create a 3-d array (replications, intervals, 2 endpoints
CIarray <- array(dim = c(replications, dim(CItemplate)),
                 dimnames = c(list(NULL), dimnames(CItemplate)))

for(n in nn) {
  cat("\nn =", n, "  ")
  set.seed(0)
  save.xbar <- rep(NA, replications)
  for(i in 1:replications) {
    if(i %in% (1:10 * (replications/10))) cat(i, " ")
    x <- rexp(n)
    v <- x / mean(x) # Compute intervals based on v, not x.
    CIarray[i, , ] <- ComputeIntervals(v, xbar = 1, n)
    save.xbar[i] <- mean(x)
  }
  # Compute E(endpoint < mu)
  gammaCI[[as.character(n)]] <- CIarray
  gammaIResults[[as.character(n)]] <- (CIarray * save.xbar > mu)
  gammaResults[[as.character(n)]] <- 1-FG(mu / pmax(CIarray, 0), n)
  # Need pmax because some endpoints are negative
}
lapply(gammaIResults, colMeans)
lapply(gammaResults, colMeans)

# Use this for standard errors
lapply(gammaResults, colStdevs)


PlotResults(gammaIResults, ylim = c(0, 0.10))
PlotResults(gammaResults, ylim = c(0, 0.10))

if(FALSE){
 save.gammaResults9 <- gammaResults
 save.gammaCI9 <- gammaCI
 save(list = c("replications", "R", "save.gammaResults9", "save.gammaCI9"),
      file = "gammaResults9.RData")
}
# save.*4 has R=10^4, 10^4 replications, and nn1
# save.*5 has R=10^4, 10^4 replications, and nn2
# save.*6 has R=10^4, 10^4 replications, nn5, and wider quantiles
# save.*7 has R=10^4, 1e3 replications, nn=c(6K, 8K), & wider quantiles
# save.*8 has R=10^4, 5e2 replications, nn=c(10K, 15K, 20K), wider quantiles
# save.*9 has R=10^4, 5e3 replications, nn=c(6K, 8K), & wider quantiles

PlotResults(c(save.gammaResults4, save.gammaResults5), ylim = c(0, 0.10))
PlotResults(c(save.gammaResults4, save.gammaResults5)[as.character(nn3)],
            ylim = c(0, 0.10))
PlotResults(c(save.gammaResults6, save.gammaResults7), ylim = c(0, 0.10))

# I'll do separate pictures instead, for better resolution.
# #pdf("../figures/coverage1.pdf", height = 4)
par(mfrow = c(1, 2), mex = .8)
PlotResults(save.normalResults6, average = TRUE, ylim = c(.02, .05),
            main = "Normal population", legend.cex = .7)
PlotResults(c(save.gammaResults4, save.gammaResults5)[as.character(nn3)],
            ylim = c(0, 0.10), main = "Exponential population",
            legend = FALSE, ylab = "")
# dev.off()




if(FALSE){
  load("normalResults4.RData")          # save.normalResults4
  load("normalResults5.RData")          # save.normalResults5
  load("normalResults6.RData")          # save.normalResults6
  load("gammaResults4.RData")           # save.gammaResults4
  load("gammaResults5.RData")           # save.gammaResults5
  load("gammaResults6.RData")           # save.gammaResults6
  load("gammaResults7.RData")           # save.gammaResults7
  load("gammaResults8.RData")           # save.gammaResults8
  load("gammaResults9.RData")           # save.gammaResults9
  for(i in seq_along(save.normalResults4))
    dimnames(save.normalResults4[[i]])[[2]][6] <- "expanded"
  for(i in seq_along(save.gammaResults4))
    dimnames(save.gammaResults4[[i]])[[2]][6] <- "expanded"
  for(i in seq_along(save.gammaResults5))
    dimnames(save.gammaResults5[[i]])[[2]][6] <- "expanded"
}

# For the arXiv version I used figures coverage{2,3,4}
# For TAS1, use coverage{2b,3b} (exclude tSkew)

# Standard error for n=8000
temp <- colStdevs(save.gammaResults9[["8000"]]) / sqrt(5000)
max(temp)
# 3.008495e-05

################
# pdf("../figures/coverage2.pdf", height = 7)
PlotResults(save.normalResults6, average = TRUE, ylim = c(.02, .05),
            main = "Normal population", legend.cex = .7, nMax = 40)
text(-1/sqrt(22), .0365, adj = 0, "p and r")
# dev.off()

# pdf("../figures/coverage3.pdf", height = 7)
PlotResults(c(save.gammaResults6, save.gammaResults9), nMax = 100,
            ylim = c(0, 0.10), main = "Exponential population")
# dev.off()
################

################ For paper2, omit tSkew (and reverse, for normal)
# pdf("../figures/coverage2b.pdf", height = 7)
PlotResults(save.normalResults6, average = TRUE, ylim = c(.02, .05),
            main = "Normal population", legend.cex = .7, nMax = 40,
            omit = c("tSkew", "reverse"))
text(-1/sqrt(22), .0365, adj = 0, "p and r")
# dev.off()

# pdf("../figures/coverage3b.pdf", height = 7)
PlotResults(c(save.gammaResults6, save.gammaResults9), nMax = 100,
            ylim = c(0, 0.10), main = "Exponential population",
            omit = "tSkew")
# dev.off()
################

################ For paper2, omit tSkew & expanded (and reverse, for normal)
# pdf("../figures/coverage2d.pdf", height = 5)
PlotResults(save.normalResults6, average = TRUE, ylim = c(.02, .05),
            main = "Normal population", legend.cex = .7, nMax = 40,
            scale = "inverse4Root",
            omit = c("tSkew", "reverse", "expanded"))
text(-1/sqrt(sqrt(22)), .0365, adj = 0, "p and r")
# dev.off()

# pdf("../figures/coverage3d.pdf", height = 7)
PlotResults(c(save.gammaResults6, save.gammaResults9), nMax = 100,
            ylim = c(0, 0.10), main = "Exponential population",
            scale = "inverse4Root",
            omit = c("tSkew", "expanded"))
# dev.off()
################

################ For paper2, combine plots
#                      omit tSkew & expanded (and reverse, for normal)
# pdf("../figures/coverage2e.pdf", height = 8)

par(mex = .8)
par(fig = c(0, 1, .62, 1), mar = c(3.1, 4.1, 4.1, 2.1))
PlotResults(save.normalResults6, average = TRUE, ylim = c(.02, .05),
            main = "Normal population", legend = TRUE, nMax = 100,
#            main = "Normal population", legend = FALSE, nMax = 100,
            scale = "inverse4Root",
            omit = c("tSkew", "reverse", "expanded"), ltyVec = c(1,4,2,3))
text(-1/sqrt(sqrt(22)), .0366, adj = 0, "(p and r)")

par(mar = c(5.1, 4.1, 4.1, 2.1))
par(fig = c(0, 1, 0, .62), new = TRUE)
PlotResults(c(save.gammaResults6, save.gammaResults9), nMax = 100,
            ylim = c(0, 0.10), main = "Exponential population",
            scale = "inverse4Root",
            omit = c("tSkew", "expanded"), ltyVec = c(1,4,2,5,3))
# dev.off()
################

################ For paper3, omit reverse
# pdf("../figures/coverage2c.pdf", height = 7)
PlotResults(save.normalResults6, average = TRUE, ylim = c(.02, .05),
            main = "Normal population", legend.cex = .7, nMax = 40,
            omit = "reverse")
text(-1/sqrt(22), .0365, adj = 0, "p and r")
# dev.off()

# pdf("../figures/coverage3c.pdf", height = 7)
PlotResults(c(save.gammaResults6, save.gammaResults9), nMax = 100,
            ylim = c(0, 0.10), main = "Exponential population",
            omit = "reverse")
# dev.off()
################


# Possible replacements for Figure 23 (showing n=5000)
combined.gammaResults <-
#  c(save.gammaResults9, save.gammaResults7, save.gammaResults8)
  c(save.gammaResults6, save.gammaResults9)
# Don't use save.gammaResults8 until I repeat that with a larger replications

PlotResults(combined.gammaResults,
            nLim = c(40, NA), scale = "inverse", nMax = 100,
            ylim = c(.01, .065), main = "Exponential population")
abline(h = c(.0225, .0275), col = "gray")
# This nicely shows the
# difference between first-order and second-order accurate

PlotResults(combined.gammaResults,
            nLim = c(40, NA), nMax = 100,
            ylim = c(.01, .065), main = "Exponential population")
abline(h = c(.0225, .0275), col = "gray")
# bootT is accurate about n=100, and tSkew about n=220. Rest after 2000.


PlotResults(combined.gammaResults,
            nLim = c(2000, NA), legend = FALSE,
            ylim = c(.021, .029), main = "Exponential population")
abline(h = c(.0225, .0275), col = "gray")
# There is more randomness apparent in the larger numbers


################
# pdf("../figures/coverage4.pdf", height = 4)
par(mfrow = c(1, 2), mex = .8)

PlotResults(combined.gammaResults,
            nLim = c(40, NA), scale = "inverse", nMax = 100, legend = FALSE,
            ylim = c(.01, .055), main = "Exponential population",
            axis2 = FALSE)
axis(side = 2, at = c(.025, .05))
abline(h = c(.0225, .0275), col = "gray")

PlotResults(combined.gammaResults,
            nLim = c(2000, NA), legend = FALSE,
            ylim = c(.021, .029), main = "Exponential population",
            axis2 = FALSE)
axis(side = 2, at = c(.0225, .025, .0275))
abline(h = c(.0225, .0275), col = "gray")
# dev.off()
################

################ for paper 3, omit reverse
# pdf("../figures/coverage4c.pdf", height = 4)
par(mfrow = c(1, 2), mex = .8)

PlotResults(combined.gammaResults,
            nLim = c(40, NA), scale = "inverse", nMax = 100, legend = FALSE,
            ylim = c(.01, .055), main = "Exponential population",
            omit = "reverse",
            axis2 = FALSE)
axis(side = 2, at = c(.025, .05))
abline(h = c(.0225, .0275), col = "gray")

PlotResults(combined.gammaResults,
            nLim = c(2000, NA), legend = FALSE,
            ylim = c(.021, .029), main = "Exponential population",
            omit = "reverse",
            axis2 = FALSE)
axis(side = 2, at = c(.0225, .025, .0275))
abline(h = c(.0225, .0275), col = "gray")
# dev.off()
################


#--------------------------------------------------
HowBigN <- function(results){
  # How big n to be 10% off for upper tail?
  means <- lapply(results, colMeans)
  means2 <- 1 - do.call("rbind", lapply(means, function(x) x[, 2]))
  n <- as.integer(rownames(means2))
  approxY <- function(...) approx(...)$y
  rbind(linear = apply(means2, 2, approxY, y = n, xout = 0.0275),
        inverse = 1/apply(means2, 2, approxY, y = 1/n, xout = 0.0275),
        inverseRt = 1/apply(means2, 2, approxY, y = 1/sqrt(n), xout = 0.0275)^2)
}
round(HowBigN(combined.gammaResults))
#              t tSkew tBoot perc reverse expanded bootT
# linear    5016   232  5260 2642      NA     2421   102
# inverse   4815   217  5063 2383      NA     2235   101
# inverseRt 4865   220  5114 2437      NA     2271   101
# Use inverse for tSkew and bootT, inverseRt for others

# When are the percentile intervals better than the t interval?
# Visually, just short of 40
HowBigNForTvsPercentile <- function(results){
  means <- lapply(results, colMeans)
  means2 <- 1 - do.call("rbind", lapply(means, function(x) x[, 2]))
  n <- as.integer(rownames(means2))
  nn <- c(20, 40, 100)
  means2 <- means2[as.character(nn), ]
  y1 <- means2[, "t"] - means2[, "perc"]
  1 / approx(y1, 1/sqrt(nn), 0)$y^2
}
HowBigNForTvsPercentile(combined.gammaResults)
# 34.08 for percentile vs t

HowBigNForTvsExpanded <- function(results){
  means <- lapply(results, colMeans)
  means2 <- 1 - do.call("rbind", lapply(means, function(x) x[, 2]))
  n <- as.integer(rownames(means2))
  nn <- c(5, 10)
  means2 <- means2[as.character(nn), ]
  y1 <- means2[, "t"] - means2[, "expanded"]
  1 / approx(y1, 1/sqrt(nn), 0)$y^2
}
HowBigNForTvsExpanded(combined.gammaResults)
# 6.36

#--------------------------------------------------
# Estimate the effective number of replications, with the variance reduction
temp.p <- do.call("rbind",
                  c(lapply(save.gammaResults4, colMeans),
                    lapply(save.gammaResults5, colMeans)))
temp.var <- do.call("rbind",
                  c(lapply(save.gammaResults4, colVars),
                    lapply(save.gammaResults5, colVars)))

temp.p <- do.call("rbind", lapply(combined.gammaResults, colMeans))
temp.var <- do.call("rbind", lapply(combined.gammaResults, colVars))
temp.replications <- sapply(combined.gammaResults, nrow)
temp.replications

# Multiplier
temp.p * (1-temp.p) / temp.var
# 8.9 for n=5, over 1000 for n=1000, over 5000 for n=8000

# Effective number of replications
temp.p * (1-temp.p) / (temp.var / rep(temp.replications, each = 7))
# at least 89K replications

#--------------------------------------------------
# sample size calculations, for howMany section

# $1.96 \sqrt{0.025\cdot 0.975/r} \le 0.025/10$,
# 1.96 sqrt(.025*.075/r) < .025/10
# 1.96^2 *.025*.975/r < .025^2/100
# r > 1.96^2 *.025*.975/(.025^2/100)
# r > 1.96^2 *.975/.025 *100
1.96^2 *.975/.025 *100
# 14982.24

# $z_{0.0275}/\zSub < s_B/\sigma_B < z_{0.0225}/\zSub$
# s_B/\sigma_B ~ N(1, var = 1/(2r)
# lhs - 1 < N(0, 1/2r) < rhs - 1  with 95% probability
# lhs - 1 < +- 1.96/sqrt(2r) < rhs - 1

qnorm(.025 * c(1.1, .9))/qnorm(.025) - 1
# -0.02096353  0.02280168
# $1.96 / \sqrt{2r} < |z_{0.0275}/\zSub - 1|$
# 1.96 / sqrt(2r) < rhs
# 1.96 / rhs < sqrt(2r)
# (1.96 / rhs)^2 / 2 < r
(1.96 / abs(qnorm(.0275)/qnorm(.025) - 1))^2 / 2
# 4370.725
(1.96 / abs(qnorm(.0225)/qnorm(.025) - 1))^2 / 2
# 3694

# Check that with simulation. Assume theoretical bootstrap population
# is normal.
x <- matrix(rnorm(10000 * 4371), 4371)
temp <- colStdevs(x)
mean(temp > qnorm(.0225)/qnorm(.025))
# .0202
mean(temp < qnorm(.0275)/qnorm(.025))
# .026
