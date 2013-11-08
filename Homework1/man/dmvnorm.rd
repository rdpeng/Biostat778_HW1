\name{dmvnorm}
\alias{dmvnorm}
\title{Fast dmv} 
\description{
 Fast calculating multivariate normal density 
}
\usage{
dmvnorm(x, mu, S, log = TRUE)
}
\arguments{
\item{x}{points to be evaluated}
\item{mu}{Normal means}
\item{S}{covariance matrix}
\item{log}{Whether log of results should be taken}
}
\details{
Make sure that S is a symmetric matrix
}
\value{
Value of normal density
}

