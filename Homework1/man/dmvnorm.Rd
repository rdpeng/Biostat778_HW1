\name{dmvnorm}
\alias{dmvnorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Evaluate multivariate normal density
%%  ~~function to do ... ~~
}
\description{This function evaluates the k-dimensional multivariate Normal density with mean mu and covariance S.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dmvnorm(x, mu, S, log = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a matrix containing the observations of multiple variable
%%     ~~Describe \code{x} here~~
}
  \item{mu}{a vector of the mean of each variable
%%     ~~Describe \code{mu} here~~
}
  \item{S}{a matrix of covariance
%%     ~~Describe \code{S} here~~
}
  \item{log}{deciding whether the output is in log form or original form
%%     ~~Describe \code{log} here~~
}
}
\details{x is a nxk matrix of points to be evaluated,mu is a length k vector of means for the k-dimensional Normal, and S is a kxk covariance matrix.
%%  ~~ If necessary, more details than the description above ~~
}
\value{The returned value is the density of these variables
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{n <- 100
n2 <- n^2
xg <- seq(0, 1, length = n)
yg <- xg
g <- data.matrix(expand.grid(xg, yg))
D <- as.matrix(dist(g))
phi <- 5

S <- exp(-phi * D)
mu <- rep(0, n2)
set.seed(1)
x <- matrix(rnorm(n2), byrow = TRUE, ncol = n2)
dmvnorm(x, mu, S, log = TRUE)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, mu, S, log = TRUE) 
{
    R = tryCatch({
        chol(S)
    }, error = function(e) {
        message("S is not positive definite")
    })
    logs = 2 * sum(log(diag(R)))
    y = x - rep(1, nrow(x)) \%*\% t(mu)
    z = forwardsolve(t(R), t(y))
    cov = diag(crossprod(z))
    final = (-length(mu)/2) * log(2 * pi) - (1/2) * logs - (1/2) * 
        cov
    if (log == TRUE) {
        final = final
    }
    else {
        final = exp(final)
    }
    return(final)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
