\name{fastlm}
\alias{fastlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit linear models fast.
}
\description{
	fastlm is used to fit linear models fast. It can be used to calculate regression coefficients, and the covariance.
}
\usage{
fastlm(x, y, na.rm = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  	a n*p matrix
}
  \item{y}{
  	a vector of length n
}
  \item{na.rm}{
  	indicated whether missing values in x or y should be removed.
}
}
\details{
This function uses the Choleski Decomposition to faster calculate the linear regression coefficients and the covariance.
}
\value{
  \item{coefficients }{a vector of the regression coefficients estimated using maximum likelihood.}
  \item{vcov }{the p*p covariance of matrix of the estimated regression coefficients.}
}

\author{
Yuan He
}

\seealso{
chol(), forwardsolve()
}
\examples{
set.seed(2)
## Generate predictor matrix
n <- 100000
p <- 500
X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))

## Coefficents
b <- rnorm(p)

## Response
y <- X %*% b + rnorm(n)

fit <- fastlm(X, y)
str(fit)
## The function is currently defined as
function (x, y, na.rm = FALSE) 
{
    if (na.rm == TRUE) {
        r = cbind(x, y)
        x = x[complete.cases(r), ]
        y = as.matrix(y[complete.cases(r)])
    }
    A = crossprod(x)
    B = crossprod(x, y)
    l = chol(A)
    beta = backsolve(l, forwardsolve(t(l), B))
    k = ncol(x)
    n = nrow(y)
    sigmaS = (crossprod(y, y) - crossprod(beta, B))/(n - k)
    covBeta = as.numeric(sigmaS) * chol2inv(l)
    list(beta, covBeta)
  }
}

\keyword{linear regression}
\keyword{covariance}
