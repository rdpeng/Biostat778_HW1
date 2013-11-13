\name{fastlm}
\alias{fastlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fast linear regression}
\description{Efficiently calculate the beta coefficients and their variances in a standard linear model: \eqn{Y = X \beta + \epsilon}}
\usage{
fastlm(X, y, na.rm = FALSE)
}
\arguments{
  \item{X}{
  The design matrix of the linear model
}
  \item{y}{
  The vector of observations
}
  \item{na.rm}{
  Should missing values be removed? If \code{TRUE}, missing values are removed. If \code{FALSE}, missing values are not removed.
}
}
\details{
This function is designed to be faster than \code{lm.fit}. It does not perform the entire suite of tasks that \code{lm} does.
}
\value{
\item{coefficients}{A vector of the estimated beta coefficients}
\item{vcov}{The variance-covariance matrix of the estimated beta coefficients}
}
\references{
Roger Peng. http://rdpeng.github.io/Biostat778_HW1/
}
\author{
Leslie Myint
}
\seealso{
\code{\link{lm}},\code{\link{lm.fit}}
}
\examples{
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
}