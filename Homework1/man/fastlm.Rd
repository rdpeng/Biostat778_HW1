\name{fastlm}
\alias{fastlm}
\title{Do fast linear models}
\usage{
  fastlm(X, y, na.rm = FALSE)
}
\arguments{
  \item{X}{a nxp matrix, where k is the number of
  measurements per observation, and n is the number of
  observations}

  \item{y}{a n-length vector of outcomes}
}
\value{
  A vector of densities, length-n
}
\description{
  \code{fastlm} Does lm, but faster. Uses cholesky
  decomposition of crossprod(X) to invert X'X.
}
\examples{
set.seed(2)
 #####Generate predictor matrix
n <- 1000
p <- 5
X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))

## Coefficents
b <- rnorm(p)

## Response
y <- X \%*\% b + rnorm(n)




system.time( fit_fast <- fastlm(X, y) )
system.time( fit_compare <- lm.fit(X, y) )
str(fit_fast)
}

