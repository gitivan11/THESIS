# Resampling regression examples

# Artificial example, bad fit
# resampling residuals
# Confidence Interval vs Prediction Interval

library(resample)

#--------------------------------------------------
# Bootstrapping a bad fit in regression (and good fit) (second version)
# This version treats the good and bad fits independently, and shows
# more variation in the good fit.

# For good fit, use skewed residuals

n <- 100
set.seed(0)
x <- runif(n, -1, 1)
z <- rnorm(n)/10
z2 <- rexp(n)/10 - .1
y <- x^2 + z
dataxy <- data.frame(x, y, z, z2)
fit <- lm(y~x)
fitz <- lm(z~x)
fitz2 <- lm(z2~x)

# Left: good fit; right: bad fit
# pdf("../figures/bootstrapGoodBadFit2.pdf", height = 3)
par(mfrow = c(1, 2), mex = .8, cex = .8, mar = rep(.1, 4))
#
# Left, Old version, normal residuals
# plot(x, z, axes = FALSE, xlab = "", ylab = "")
# bootstrap(dataxy, abline(lm(z~x, data = dataxy)), R = 20, seed = 1)
# abline(fitz, col = "white", lwd = 3)
# abline(fitz, lwd = 3, lty = 2, col = "red")
# Left, New version, skewed residuals
plot(x, z2, axes = FALSE, xlab = "", ylab = "")
bootstrap(dataxy, abline(lm(z2~x, data = dataxy)), R = 20, seed = 1)
abline(fitz2, col = "white", lwd = 3)
abline(fitz2, lwd = 3, lty = 2, col = "red")
# right: bad fit
plot(x, y, ylim = range(y,z), axes = FALSE, xlab = "", ylab = "")
bootstrap(dataxy, abline(lm(y~x, data = dataxy)), R = 20, seed = 1)
abline(fit, col = "white", lwd = 3)
abline(fit, lwd = 3, lty = 2, col = "red")
legend("top", lwd = 2:1, col = c("red", "black"), lty = 2:1,
       legend = c("Original line", "Bootstrap lines"))
# dev.off()


#--------------------------------------------------
# (1) resampling residuals

set.seed(0)
x1 <- 1:20 - runif(20)
y1 <- x1 + rnorm(20) * 4
plot(x1, y1, axes = FALSE, xlab = "", ylab = "")
fit1 <- lm(y1 ~ x1)
abline(fit1)
segments(x1, y1 = y1, y0 = predict(fit1))
