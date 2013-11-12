\name{dmvnorm}
\alias{dmvnorm}
\title{get multivariate normal density for a matrix of vectors.}
\usage{
  dmvnorm(x, mu, S, log = TRUE)
}
\arguments{
  \item{x}{nxk matrix, where k is the number of
  measurements per observation, and n is the number of
  observations}

  \item{mu}{k-length vector of means}

  \item{S}{a kxk covariance matrix}

  \item{log}{whether to output log denisites or not.}
}
\value{
  A vector of densities, length-n
}
\description{
  \code{dmvnorm} gets the multivariate density.
}
\examples{
kSqrt <- 5 #our final density will be for a random vector of length 15
k <- kSqrt^2
xg <- seq(0, 1, length = kSqrt)
yg <- xg
g <- data.matrix(expand.grid(xg, yg))
D <- as.matrix(dist(g))
phi <- 5
S <- exp(-phi * D) #S ends up being k by k
n<-4 #four observations

mu <- rep(0, k)
set.seed(1)
x <- matrix(rnorm(k*n), byrow = TRUE, ncol = k) #n x k, this is the usual input to dmvnorm i mvtnorm

system.time( fastOutput <-dmvnorm(x=x,mu=mu,S=S,log=TRUE) )
str(fastOutput)
}

