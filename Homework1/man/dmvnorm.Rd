\name{dmvnorm}
\alias{dmvnorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast Multivariate Normal Density
}
\description{
This function evaluates the k-dimensional multivariate Normal density with mean mu and covariance S.
}
\usage{
dmvnorm(x, mu, S, log = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
 a n*k matrix of points to be evaluated.
}
  \item{mu}{
 a vector of means of length k for the k-dimensional Normal.
}
  \item{S}{
 a k*k covariance matrix.
}
  \item{log}{
 If log == TRUE, returns the logged density(by default), otherwise, returns the original density.
}
}
\details{
 In this function, Cholesky decomposition was used to compute the Mahalanobis distance term and the determinant of S in the density function. By using the function chol(), the positive definiteness of S is checked. If it is not, the function will return a error message.
}
\value{
The function returns a vector of length n, containing the values of the multivariate Normal density evaluated at the n points.
}
\references{
Nocedal, Jorge, and S. Wright. Numerical optimization, series in operations research and financial engineering. Springer, New York (2006).
}
\author{
Detian Deng
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{
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

mymvpdf<-dmvnorm(x=x, mu=mu, S=S, log = TRUE)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line