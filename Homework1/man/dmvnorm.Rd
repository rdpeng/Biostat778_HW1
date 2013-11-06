\name{dmvnorm}
\alias{dmvnorm}
\title{Evaluate the density of a multivariate normal distribution}
\usage{
  dmvnorm(x, mu, S, log = TRUE)
}
\arguments{
  \item{x}{an \code{n x k} matrix to represent \code{n}
  vectors of length \code{k}}

  \item{mu}{the mean of the distribution: a numeric vector
  of length \code{k}}

  \item{S}{the variance-covariance matrix with dimensions
  \code{k x k}}

  \item{log}{a boolean that indicates if the log of the
  density should be returned}
}
\value{
  numeric vector with length \code{n}, one entry for the
  density of the \code{n} vectors in \code{x} evaluated
  under a multivariate normal distribution defined by mean
  \code{mu} and variance-covariance \code{S}
}
\description{
  Evaluate the density of a multivariate normal
  distribution
}
\examples{
x <- matrix(rnorm(10*9), ncol=9)
mu <- rep(0,9)
xg <- seq(0, 1, length = 3)
yg <- xg
g <- data.matrix(expand.grid(xg, yg))
D <- as.matrix(dist(g))
S <- exp(D * -1)
dmvnorm(x, mu, S)
}

