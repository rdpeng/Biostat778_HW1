\name{fastlm}
\alias{fastlm}
\title{Fast Linear Model Fit}
\description{
This function fast fits a linear regression model using cholesky decomposition. Fitted beta and its covariance are returned.
}
\usage{
fastlm(X, y, na.rm=FALSE)
}
\arguments{
\item{X}{a n by p matrix.}
\item{y}{a vector of length n.}
\item{na.rm}{a logical variable that indicates whether missing values in X or y should be removed.}

}
\keyword{linear model,multivariate normal,fast}
\examples{
set.seed(2)
## Generate predictor matrix
n <- 1000
p <- 500
X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))

## Coefficents
b <- rnorm(p)

## Response
y <- X %*% b + rnorm(n)

fit <- fastlm(X, y)
str(fit)
}