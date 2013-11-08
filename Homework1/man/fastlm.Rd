\name{fastlm}
\alias{fastlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Fit a linear regression model
}
\description{
  Fit a linear regression model using fast computation
}
\usage{
fastlm(X, y, na.rm = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
    The n x p matrix of p predictors for n observations
}
  \item{y}{
    The length n vector of responses
}
  \item{na.rm}{
    Should missing values from \code{X} or \code{y} be removed first?
    Default is \code{FALSE}.
}
}
\details{
  This function computes the linear regression coefficients and the covariance
  matrix of the estimated coefficients in a fast manner.
}
\value{
A list containing the following items is returned:
 \item{coefficients}{A length p vector of estimated regression coefficients}
 \item{vcov}{A p x p matrix containing the covariance matrix of the estimated
   regression coefficients}
}
\author{
  Roger D. Peng \email{rpeng@jhsph.edu}
}
\keyword{models}
