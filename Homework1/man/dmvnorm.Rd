\name{dmvnorm}
\alias{dmvnorm}
\title{Fast Multivariate Normal (log)density Calculation }
\description{
This function fast calculates the density(or the log density) of a multivariate normal distribution, using cholesky decomposition. A vector of length n is returned.
}
\usage{
dmvnorm(x, mu, S, log = TRUE)
}
\arguments{
\item{x}{a n by k matrix of points to be evaluated.}
\item{mu}{a length k vector; mean of mvn distribution.}
\item{S}{a k by k covariance matrix.}
\item{log}{a logical variable that indicates whether the log density should be returned.}
}

\keyword{multivariate normal,density,fast}
\examples{
## Create the covariance matrix
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

dmvnorm(x, mu, S, log = TRUE)
}