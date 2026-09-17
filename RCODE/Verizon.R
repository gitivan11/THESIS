# Analyze Verizon data

library(resample)
library(Hmisc) # For latex function

#--------------------------------------------------

data(Verizon) # This comes with library(resample)
# Verizon <- read.csv("../data/Verizon.csv")

ILEC <- with(Verizon, Time[Group == "ILEC"])
CLEC <- with(Verizon, Time[Group == "CLEC"])
cbind(ILEC = c(n = length(ILEC), mean = mean(ILEC), sd = sd(ILEC)),
      CLEC = c(length(CLEC), mean(CLEC), sd(CLEC)))
#             ILEC     CLEC
# n    1664.000000 23.00000
# mean    8.411611 16.50913
# sd     14.690039 19.50358
mean(CLEC) - mean(ILEC)
# 8.09752


#system("mkdir ../figures")
# pdf("../figures/VerizonData.pdf", height = 4, width = 4)
par(mex = .6, mar = c(5.1, 5.1, 2.1, .1))
qqnorm(ILEC, ylab = "Repair Times (hours)", main = "")
points(qqnorm(CLEC, plot.it = FALSE), col = 2, pch = 2)
legend(x = "topleft", pch = 1:2, col = 1:2, legend = c("ILEC", "CLEC"))
# dev.off()


#### Permutation test
permVerizon <-
    permutationTest2(Verizon, mean(Time), treatment = Group,
                     alternative = "greater", seed = 0,
                     statisticNames = "mean")
permVerizon
# P-value 0.0171

hist(permVerizon, xlab = "Difference in means, CLEC-ILEC (hours)",
     ylab = "Verizon data")

#### Two-sample bootstrap test by pooling.
# The resample package does not support this.
set.seed(0)
r <- 10^4 - 1
replicates1 <- replicates2 <- numeric(r)
i1 <- samp.bootstrap(length(Verizon$Time), r, size = length(CLEC))
i2 <- samp.bootstrap(length(Verizon$Time), r, size = length(ILEC))
for(i in 1:r) {
  replicates1[i] <- mean(Verizon$Time[i1[, i]])
  replicates2[i] <- mean(Verizon$Time[i2[, i]])
}
replicatesDiff <- replicates1 - replicates2
(1+sum(replicatesDiff >= permVerizon$observed)) / (r+1)
# P-value = .0183. Slightly larger than perm, within random variation
hist(replicatesDiff)
hist.resample(list(p = 1, replicates = matrix(replicatesDiff),
                   observed = mean(CLEC) - mean(ILEC),
                   stats = list(Mean = mean(replicatesDiff))),
              xlab = "Difference in means, CLEC-ILEC (hours)",
              ylab = "Pooled Bootstrap Test")


#### t-tests
# P-values from t tests
t.test(Time ~ Group, data = Verizon, alternative = "greater")
# t = 1.9834, df = 22.346, p-value = 0.02987
t.test(Time ~ Group, data = Verizon, alternative = "greater", var.equal = TRUE)
# t = 2.6125, df = 1685, p-value = 0.004534

# In real life, the Public Utilities Commision mandated the use
# of a t-statistic with denominator calculated solely from the ILEC group.
# This makes the standard error artificially low, when the CLEC group
# has larger mean (and larger SE), and makes the P-value even smaller
t1 <- (mean(CLEC) - mean(ILEC)) /
    (sd(ILEC) * sqrt(1/length(ILEC) + 1/length(CLEC)))
df1 <- length(ILEC) - 1
cat(sep = "", "t = ", round(t1, 4), ", df = ", round(df1, 3),
    ", p-value = ", round(pt(-t1, df = df1), 5), "\n")
# t = 2.6255, df = 1663, p-value = 0.00437


#### Bootstrap
bootCLEC <- bootstrap(CLEC, mean, seed = 0)
bootILEC <- bootstrap(ILEC, mean, seed = 1)
bootCLEC
#       Observed       SE     Mean      Bias
# stat1 16.50913 3.961816 16.53088 0.0217463
bootILEC
#       Observed       SE     Mean         Bias
# stat1 8.411611 0.357599 8.404107 -0.007503157


CI.percentile(bootILEC, expand = FALSE)
#          2.5%    97.5%
# mean 7.728359 9.129599

CI.percentile(bootCLEC, expand = FALSE)
#          2.5%    97.5%
# mean 10.09124 25.40154

CI.percentile(bootCLEC, expand = FALSE) - bootCLEC$observed
# mean -6.417891 8.892413

2 * bootCLEC$observed - CI.percentile(bootCLEC, expand = FALSE) # reverse percentile
# mean 22.92702 7.616717

t.test(ILEC)
t.test(CLEC)

# pdf("../figures/bootVerizon1.pdf")
par(mfrow = c(2, 2), mex = .8)
hist(bootILEC, xlab = "Commercial Time (ILEC)")
qqnorm(bootILEC, pch = ".")
qqline(bootILEC$replicates)
hist(bootCLEC, xlab = "Commercial Time (CLEC)")
qqnorm(bootCLEC, pch = ".")
qqline(bootCLEC$replicates)
# dev.off()
par(mfrow = c(1, 1), mex = 1)

# Monte Carlo SE for the endpoints of the bootstrap percentile interval
CI.percentile(bootCLEC, expand = FALSE)
# mean 10.09124 25.40154
bootMC <- bootstrap(bootCLEC$replicates,
                    quantile(data, probs = c(.025, .975), type = 6), seed = 2)
bootMC
#       Observed         SE     Mean        Bias
# 2.5%  10.09124 0.06628435 10.09870 0.007458429
# 97.5% 25.40154 0.14100334 25.41162 0.010081198

CI.percentile(bootCLEC, expand = FALSE) # matches Observed

#### Two-sample bootstrap
bootVerizon2 <- bootstrap2(Verizon, mean(Time), treatment = Group,
                           seed = 0,
                           statisticNames = "Difference in means")
bootVerizon2
#                     Observed       SE     Mean       Bias
# Difference in means  8.09752 3.979942 8.113467 0.01594722


CI.percentile(bootVerizon2, expand = FALSE)
#                         2.5%    97.5%
# Difference in means 1.632452 16.99053
CI.percentile(bootVerizon2, probs = c(.01), expand = FALSE)
#                                       1%
# Difference in means: CLEC-ILEC 0.7718843

# pdf("../figures/bootVerizon3.pdf", height = 4)
par(mfrow = c(1, 2), mex = .8)
hist(bootILEC, xlab = "Commercial Time (ILEC)")
hist(bootCLEC, xlab = "Commercial Time (CLEC)")
# dev.off()
par(mfrow = c(1, 1), mex = 1)

# Two-sample bootstrap and permutation distributions.
par(mfrow = c(1, 2), mex = .8)
hist(bootVerizon2,
     ylab = "Verizon data")
hist(permVerizon, xlab = "Difference in means, CLEC-ILEC (hours)",
     ylab = "Verizon data")
par(mfrow = c(1, 1), mex = 1)

# pdf("../figures/bootPermVerizon2.pdf")
par(mfrow = c(2, 2), mex = .8)
hist(bootVerizon2, xlab = "Bootstrap (CLEC-ILEC)")
qqnorm(bootVerizon2, pch = ".")
qqline(bootVerizon2$replicates)
hist(permVerizon, xlab = "Permutation Test (CLEC-ILEC)")
qqnorm(permVerizon, pch = ".")
qqline(permVerizon$replicates)
# dev.off()
par(mfrow = c(1, 1), mex = 1)


#### Bootstrap t for ILEC
bootT.ILEC <- bootstrap(ILEC, c(mean(ILEC), sd(ILEC)), seed = 1)
ILECt <- with(bootT.ILEC,
              (replicates[, 1]-observed[1])/(replicates[, 2]/sqrt(n)))

c(mean(ILECt < qt(.025, 1664-1)),
  mean(ILECt > qt(.975, 1664-1)))
# 0.0339, 0.0166

# Takes a while (higher accuracy)
bootT.ILEC6 <- bootstrap(ILEC, c(mean(ILEC), sd(ILEC)), seed = 1, R = 10^6)
ILECt6 <- with(bootT.ILEC6,
               (replicates[, 1]-observed[1])/(replicates[, 2]/sqrt(n)))
c(mean(ILECt6 < qt(.025, 1664-1)),
  mean(ILECt6 > qt(.975, 1664-1)))
# 0.035878, 0.016931

# pdf("../figures/bootVerizonT.pdf", height = 4)
par(mfrow = c(1, 2), mex = .8)
# Scatterplot
plot(bootT.ILEC$replicates[, 1], bootT.ILEC$replicates[, 2]/sqrt(bootT.ILEC$n),
     xlab = "Mean", ylab = "SE", pch = ".")
# bootstrap t distribution
plot(qt(ppoints(10^4), 1664-1), sort(ILECt),
     pch = ".", xlab = "Theoretical t Quantiles",
     ylab = "Bootstrap t Quantiles", main = "T-T Quantile Plot")
abline(0, 1)
# dev.off()
par(mfrow = c(1, 1), mex = 1)



# Compare skewness of bootstrap mean and bootstrap t
par(mfrow = c(1, 2), mex = .8)
qqnorm(bootT.ILEC, resampleColumns = 1, ylab = "Bootstrap Means")
abline(mean(ILEC), sd(ILEC)/sqrt(length(ILEC)), lty=2)
plot(qt(ppoints(10^4), 1664-1), sort(ILECt),
     pch = ".", xlab = "Theoretical t Quantiles",
     ylab = "Bootstrap t Quantiles", main = "T-T Quantile Plot")
abline(0, 1, lty=2)
par(mfrow = c(1, 1), mex = 1)


# For table showing endpoints, above and below xbar
temp <- data.frame(bootstrapT = CI.bootstrapT(bootT.ILEC)[1,],
                   tWithBootSE = CI.t(bootT.ILEC)[1, ],
                   percentile = CI.percentile(bootT.ILEC, expand = FALSE)[1, ])
temp
#       bootstrapT tWithBootSE percentile
# 2.5%    7.765737    7.710219   7.728359
# 97.5%   9.173910    9.113002   9.129599
temp - bootT.ILEC$observed[1]
#       bootstrapT tWithBootSE percentile
# 2.5%  -0.6458736  -0.7013917 -0.6832520
# 97.5%  0.7622997   0.7013917  0.7179887
latex(round(temp[, c(2, 3, 1)] - bootT.ILEC$observed[1], 3), file = "")
temp2 <- temp[, c(2, 3, 1)] - bootT.ILEC$observed[1]
temp2[2,] / temp2[1,]
#       tWithBootSE percentile bootstrapT
# 97.5%          -1   -1.05084  -1.180261
.765 / .648 # tSkew, also 1.180

# What does the bootstrap t say about the actual coverage of the other two?
temp2 <- with(bootT.ILEC, (replicates[, 1] - observed[1])/replicates[, 2])
temp3 <- as.matrix(temp)
for(i in 1:6)
  temp3[i] <- with(bootT.ILEC,
                   1-mean(temp2 < (observed[1] - temp3[i])/observed[2]))
temp3
#       bootstrapT tWithBootSE percentile
# 2.5%       0.025      0.0172     0.0199
# 97.5%      0.975      0.9652     0.9681
1 - temp3[2, ]
#      0.0250      0.0348      0.0319

temp4 <- rbind(temp3[1, ], 1 - temp3[2, ])
latex(round(temp4[, c(2, 3, 1)], 3), file = "")
round(temp4[, c(2, 3, 1)] / .025, 2)
#      tWithBootSE percentile bootstrapT
# [1, ]        0.69       0.80          1
# [2, ]        1.39       1.28          1



Skewness <- function(x){
  x <- x-mean(x)
  mean(x^3) / mean(x^2)^1.5
}
# Skewness(qexp(ppoints(10^4))) # 1.99
# Skewness(qgamma(ppoints(10^4), shape = 4)) # .997
# Skewness of gamma is 2/sqrt(shape)


# Asymptotic result - Edgeworth expansion to estimate noncoverages
# of a t interval for the Verizon data.
# kappa <- skewness/(6 * sqrt(n))
# alpha2 <- c(.025, .975)
# talpha <- qt(alpha2, n-1)
# alpha2 + (2 * talpha^2 + 1) * dnorm(talpha)
nILEC <- length(ILEC)
kappa <- Skewness(ILEC)/(6 * sqrt(nILEC))
alpha2 <- c(.025, .975)
talpha <- qt(alpha2, nILEC-1)
alpha2 + kappa * (2 * talpha^2 + 1) * dnorm(talpha)
# 0.03446471 0.98446471   # 1 - .98446471 = 0.01553529
alpha2 + kappa * (2 * talpha^2 + 1) * dt(talpha, nILEC-1)
# practically the same. dt / dnorm = 1.000916, so little difference there.

# How large n to reduce errors to 10%?
# n \ge \left( \frac{\gamma}{6} \frac{10}{\alpha}
#   (2 \zSub^2 + 1) \phi(\zSub) \right)^2


(Skewness(ILEC)/6 * 10/0.025 * (2 * qnorm(.975)^2 + 1) * dnorm(qnorm(.975)))^2
# 23922.12

#--------------------------------------------------
# bias of s^2 and sigmahat^2

Sigmahat <- function(x) {
  varhat <- mean((x-mean(x))^2)
  c("s^2" = var(x),
    s = sd(x),
    varhat = varhat,
    sigmahat = sqrt(varhat))
}

bootSigmahat <- bootstrap(CLEC, Sigmahat, seed = 0)
bootSigmahat
#           Observed         SE      Mean       Bias
# s^2      380.38946 267.201641 362.36533 -18.024129
# s         19.50358   7.352757  17.55869  -1.944883
# varhat   363.85079 255.584179 346.61032 -17.240471
# sigmahat  19.07487   7.191139  17.17274  -1.902133

var(CLEC)
# 380.3895
-Sigmahat(CLEC)[2]/length(CLEC)
# -15.8196
# The Bias estimate is larger, but this is random chance; with R=10^6 is -16.04

CI.percentile(bootSigmahat, expand = FALSE)
#               2.5%     97.5%
# s^2      59.398073 931.29641
# s         7.707015  30.51715
# varhat   56.815548 890.80526
# sigmahat  7.537609  29.84636

CI.t(bootSigmahat)
#                 2.5%     97.5%
# s^2      -173.752825 934.53175
# s           4.254890  34.75226
# varhat   -166.198354 893.89994
# sigmahat    4.161365  33.98838

# T intervals go negative; could use log transform.

CI.t(bootSigmahat)[2,]^2

# Classical chi-square interval for sigma^2
# (n-1)s^2 / qchisq(.975, 22) < sigma^2 < (n-1)s^2 / qchisq(.025, 22)
22 * var(CLEC) / qchisq(c(.975, .025), 22)
# 227.5260 762.0036

# That assumes (n-1)s^2/sigma^2 ~ chisq(n-1)
#  var(s^2) = s^4/(n-1)^2 * 2(n-1) = 2 s^4/(n-1)
sqrt(2 * var(CLEC)^2 / 22)
# 114.7 - compare to bootstrap SE^2 of 267
