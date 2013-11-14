\name{dmvnorm}
\alias{dmvnorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   multivariate normal density 
}
\description{
Quickly evaluate multivariate normal density.
}
\usage{
dmvnorm(x, mu, S, log = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{mu}{
%%     ~~Describe \code{mu} here~~
}
  \item{S}{
%%     ~~Describe \code{S} here~~
}
  \item{log}{
%%     ~~Describe \code{log} here~~
}
}
\details{
    Uses cholesky transformation with backwards/forwards substitution to find
    coefficient estimates and covariance of beta_hat.
}
\author{
    Stephen Cristiano
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, mu, S, log = TRUE) 
{
    if (!is.matrix(x)) 
        x <- matrix(x, nrow = 1, ncol = ncol(S))
    k <- ncol(S)
    U <- try(chol(S), silent = TRUE)
    if (class(U) == "try-error") 
        stop("S is not positive definite")
    d <- diag(U)
    logd <- sum(log(d))
    b <- crossprod(forwardsolve(t(U), t(x - mu)))
    if (is.matrix(b)) 
        b <- diag(b)
    logf <- -k/2 * log(2 * pi) - logd - 0.5 * b
    if (log) 
        return(logf)
    else return(exp(logf))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
