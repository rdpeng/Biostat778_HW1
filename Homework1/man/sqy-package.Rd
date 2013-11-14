\name{sqy-package}
\alias{sqy-package}
\alias{sqy}
\docType{package}
\title{
What the package does (short line)
~~ package title ~~
}
\description{
The package contains 2 functions.One is for fast linear model regression and another for Multivariate normal distribution density.}
\details{
\tabular{ll}{
Package: \tab sqy\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-11-13\cr
License: \tab What license is it under?\cr
}
~~ An overview of how to use the package, including the most important functions ~~
}
\author{JerryPPA(Qiyuan Shi)
Who wrote it

Maintainer: Who to complain to <yourfault@somewhere.net>
~~ The author and/or maintainer of the package ~~
}
\references{
~~ Literature or other references for background information ~~
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation ~~
~~ directory ~~
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
%
n <- 100
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
%examples for fastlm
set.seed(2)
n <- 100000
p <- 500
X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
b <- rnorm(p)
y <- X %*% b + rnorm(n)
fit <- fastlm(X, y)
str(fit)


}
