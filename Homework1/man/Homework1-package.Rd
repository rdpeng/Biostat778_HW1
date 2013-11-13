\name{YuanHeHW1-package}
\alias{YuanHeHW1-package}
\alias{YuanHeHW1}
\docType{package}
\title{
Fast linear regression
}
\description{
This is a mini package that realise fast calculation of linear regression coefficients and multivariate normal density.
}
\details{
\tabular{ll}{
Package: \tab Yuan_HW1\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-11-11\cr
License: \tab GPL \cr
}
}
\author{
Yuan He

Maintainer: Yuan He <yhe23@jhu.edu>
~~ The author and/or maintainer of the package ~~
}
\references{
~~ Literature or other references for background information ~~
}
\keyword{fastRpackages}
\examples{
set.seed(2)
n <- 100000
p <- 500
x <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
b <- rnorm(p)
y = x %*% b + rnorm(n)

fit <- fastlm(x, y)
str(fit)


## Create the covariance matrix
n <- 3
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

dmvnorm(x, mu, S, log = TRUE)
}
