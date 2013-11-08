\name{Homework1-package}
\alias{Homework1-package}
\alias{Homework1}
\docType{package}
\title{
HW1 fastlm and dmvnorm
}
\description{
This package is for 140.778 HW1. It has two functions. 
The fastlm function fits a linear model given the design matrix X and
outcome vector Y. The dmvnorm function evaluates the multivariate normal 
density given the matrix of quantiles, vector of mean, and covariance matrix. 
}
\details{
\tabular{ll}{
Package: \tab Homework1\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-11-05\cr
License: \tab What license is it under?\cr
}
This package has two functions. 
The 'fastlm' function fits a linear model given the design matrix X and
outcome vector Y. 
The 'dmvnorm' function evaluates the multivariate normal 
density given the matrix of quantiles, vector of mean, and covariance matrix. 
}
\author{
Detian Deng

Maintainer: Detian Deng <ddeng@jhsph.edu>
}
\references{
Nocedal, Jorge, and S. Wright. "Numerical optimization, series in operations research and financial engineering." Springer, New York (2006).
}

%\keyword{ package }

%\seealso{
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}

\examples{
# fastlm
set.seed(2)
## Generate predictor matrix
n <- 100
p <- 15
X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))

## Coefficents
b <- rnorm(p)

## Response
y <- X \%*\% b + rnorm(n)

fit <- fastlm(X, y)
str(fit)

# dmvnorm
n <- 5
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

mymvpdf<-dmvnorm(x=x, mu=mu, S=S, log = TRUE)

}
