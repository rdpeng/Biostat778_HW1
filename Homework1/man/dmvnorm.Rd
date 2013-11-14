\name{dmvnorm}
\alias{dmvnorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast Multivariate Normal Density
}
\description{
This function evaluates the k-dimensional multivariate normal desnity
with mean mu and covariance S. S must be positive definite and symmetric.
The function returns a vector of length n containing the values of the 
multivariate normal density evaluated at n points. 
}
\usage{
dmvnorm(x, mu, S, log = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a nxk matrix of points to be evaluated
}
  \item{mu}{
a length k vector of means for the k-dimensional Normal
}
  \item{S}{
a kxk covariance matrix 
}
  \item{log}{
logical. If log=TRUE then the function returns the log density.
}
}
\details{
The covariance matrix S must be symmetric. The function checks if it is
positive definite. 
}
\value{
A vector of length n containing the values of the multivariate normal 
density evaluated at n points.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Shan Andrews
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
op <- options(scipen = 5)

n <- 10
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

y <- dmvnorm(x, mu, S, log = TRUE)
print(y, digits = 10)

options(op)

## The function is currently defined as
function (x, mu, S, log = TRUE) 
{
    checkem <- eigen(S)$values
    if (sum(checkem > 0) != length(checkem)) {
        stop("S is not positive definite")
    }
    k <- length(x)
    expterm <- ((-1/2) * sqrt(mahalanobis(x, mu, S)))
    dens <- exp(((-k/2) * log(2 * pi)) - ((1/2) * log(det(S))) + 
        expterm)
    if (log == TRUE) {
        return(log(dens))
    }
    return(dens)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
