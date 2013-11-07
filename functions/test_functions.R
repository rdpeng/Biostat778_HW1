setwd("~/Documents/Johns Hopkins SPH/Git Repo/Biostat778_HW1/functions")

set.seed(2)
## Generate predictor matrix
n <- 100000
p <- 500
X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))

## Coefficents
b <- rnorm(p)

## Response
y <- X %*% b + rnorm(n)

source("fastlm.R")
system.time(fit <- fastlm(X, y))
str(fit)

########################################################
n <- 50
n2 <- n^2
xg <- seq(0, 1, length = n)
yg <- xg
g <- data.matrix(expand.grid(xg, yg))
D <- as.matrix(dist(g))
phi <- 5

S <- exp(-phi * D)
mu <- rep(0, n2)
set.seed(1)
x <- matrix(rnorm(n2), byrow = TRUE, ncol = n2)

source("fastdmvnorm.R")
system.time(mymvpdf<-dmvnorm(x=x, mu=mu, S=S, log = TRUE))

library(mvtnorm)
system.time(stdmvpdf<-dmvnorm(x=x,mean=mu,sigma=S,log=TRUE))


#
package.skeleton(name="Detian_HW1",list=c("fastlm","dmvnorm"))



