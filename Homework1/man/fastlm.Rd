\name{fastlm}
\alias{fastlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast Linear Regression
}
\description{
fastlm() fits a linear regression model to outcome data y and predictor data in a matrix X.
The inputs should be X, a n*p matrix, y, a vector of length n, and na.rm, which indicates whether missing values in X or y should be removed.
This function can always run faster than the lm.fit() function in R using the same inputs.}
\usage{
fastlm(X, y, na.rm = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
a n*p matrix, i.e. the design matrix with n observations and (p-1) features
}
  \item{y}{
a vector of length n, i.e. the n obervations of outcome variable
}
  \item{na.rm}{
if na.rm == TRUE, then remove observations with NA, otherwise, keep them.
}
}
\details{
In this function, Cholesky decomposition was applied to solve the normal equation and to compute the inverse of t(X) X, which was needed for estimating covariance matrix of the estimated regression coefficients. 

}
\value{
\item{coefficients}{ a vector of the regression coefficients estimated using maximum likelihood}

\item{vcov}{ the p*p covariance matrix of the estimated regression coefficients.}

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
set.seed(2)
## Generate predictor matrix
n <- 100
p <- 15
X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))

## Coefficents
b <- rnorm(p)

## Response
y <- X \%*\% b + rnorm(n)

fit <- fastlm(X, y)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line