\name{Homework1-package}
\alias{Homework1-package}
\alias{Homework1}
\docType{package}
\title{
Fast linear regression and multivariate normal density calculation
~~ package title ~~
}
\description{
More about what it does (maybe more than one line)
~~ A concise (1-5 lines) description of the package ~~
}
\details{
\tabular{ll}{
Package: \tab Homework1\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-11-13\cr
License: \tab NA\cr
}
This package contains two functions. The first function fits a linear regression model to outcome data y and predictor data in a matrix x. The second function evaluates a k-dimensional multivariate Normal density with mean mu and covariance S. 
}
\author{
Shan Andrews

Maintainer: Shan Andrews <sandrews@jhsph.edu>
~~ The author and/or maintainer of the package ~~
}
\references{
None
}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
%~~ the R documentation directory ~~
\keyword{ regression, normal  distribution }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
op <- options(scipen = 5)

n <- 10
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

y <- dmvnorm(x, mu, S, log = TRUE)
print(y, digits = 10)

options(op)

op <- options(scipen = 5, width = 80)

set.seed(2)
## Generate predictor matrix
n <- 1000
p <- 50
X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))

## Coefficents
b <- rnorm(p)

## Response
y <- X %*% b + rnorm(n)

fit <- fastlm(X, y)
print(drop(fit$coefficients))
diag(fit$vcov)

options(op)
}
