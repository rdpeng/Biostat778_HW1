\name{dmvnorm}
\alias{dmvnorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Compute the multivariate Normal density
}
\description{
  Compute a fast version of the multivariate Normal density function
}
\usage{
dmvnorm(x, mu, S, log = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    The n x k matrix of n points at which to evaluate the density 
}
  \item{mu}{
    The length k vector of means
}
  \item{S}{
    The k x k covariance matrix
}
  \item{log}{
    Should the function return the log density? Default is \code{TRUE}.
}
}
\details{
  This function checks to make sure that \code{S} is positive definite
  and returns an error if it is not.
}
\value{
  A length n vector of (log) density values is returned
}
\author{
  Roger D. Peng \email{rpeng@jhsph.edu}
}
\keyword{distribution}
