\name{fastlm}
\alias{fastlm}
\title{Fast linear regression}
\usage{
  fastlm(X, y, na.rm = FALSE)
}
\arguments{
  \item{X}{an \code{n x p} design matrix}

  \item{y}{a numeric vector of length \code{n}}

  \item{na.rm}{a boolean to indicate if \code{NA} values
  should be removed}
}
\value{
  \code{fastlm()} returns a list with the following
  components: \item{coefficients}{A numeric vector of
  length \code{p} containing the regression coefficients
  estimated by maximum likelihood.} \item{vcov}{the \code{p
  x p} covariance matrix of the estimated regression
  coefficients.}
}
\description{
  Fast linear regression
}
\examples{
X <- matrix(rnorm(10),ncol=2)
y <- rnorm(5)
fastlm(X, y)
}

