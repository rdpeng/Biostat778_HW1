\name{dmvnorm}
\alias{dmvnorm}

\title{
Fast Multivariate Normal Density Calculation
}
\description{
Return a vector of length n containing the values of the multivariate Normal density evaluated at the n points
}
\usage{
dmvnorm(x, mu, S, log = TRUE)
}
\arguments{
  \item{x}{
an n by k matrix of points to be evaluated
}
  \item{mu}{
a length k vector of means
}
  \item{S}{
a k by k covariance matrix
}
  \item{log}{
If true,  return the log density.
}
}