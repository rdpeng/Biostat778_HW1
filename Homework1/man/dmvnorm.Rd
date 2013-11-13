\name{dmvnorm}
\alias{dmvnorm}
\title{Fast multivariate normal density}
\description{Fast evaluation of the multivariate normal density}
\usage{
dmvnorm(x, mu, S, log = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  The \eqn{n} by \eqn{k} matrix of points at which to evaluate the density. Each row of the matrix corresponds to one point (of dimension \eqn{k}).
}
  \item{mu}{
  The length \eqn{k} vector of means of the multivariate normal.
}
  \item{S}{
  The \eqn{k} by \eqn{k} covariance matrix of the multivariate normal.
}
  \item{log}{
  If \code{TRUE}, returns the natural logarithm of the density. If \code{FALSE}, returns the density.
}
}
\details{
This function checks that \code{S} is positive-definite. If not, the function will stop.
}
\value{
A vector of length \eqn{n} containing the multivariate normal evaulated at the \eqn{n} points of \code{x}.
}
\references{
Roger Peng. http://rdpeng.github.io/Biostat778_HW1/
}
\author{
Leslie Myint
}
\seealso{
\code{\link{mvrnorm}}
}
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