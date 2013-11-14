\name{fastlm}
\alias{fastlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast Linear regression
}
\description{
This function fits a linear regression model to outcome data y and
predictor data in a matrix X. 
}
\usage{
fastlm(X, y, na.rm = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
An nXp matrix of predictor data
}
  \item{y}{
A vector of length n of outcome data
}
  \item{na.rm}{
Logical indicating whethere missing values in X or y should be removed
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
	\item{coefficients }{A vector of regression coefficients estimated using maximum likelihood}
	\item{vcov }{A pxp covariance matrix of the estimated regression coefficients}
%% ...
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
op <- options(scipen = 5, width = 80)

set.seed(2)
## Generate predictor matrix
n <- 1000
p <- 50
X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))

## Coefficents
b <- rnorm(p)

## Response
y <- X %*% b + rnorm(n)

fit <- fastlm(X, y)
print(drop(fit$coefficients))
diag(fit$vcov)

options(op)

## The function is currently defined as
function (X, y, na.rm = FALSE) 
{
    if (na.rm == TRUE) {
        Xrm <- which(is.na(X), arr.ind = TRUE)
        yrm <- which(is.na(y))
        getrid <- sort(unique(c(Xrm[, 1], yrm)))
        X <- X[-getrid, ]
        y <- y[-getrid]
    }
    first <- as.matrix(solve(t(X) \%*\% X))
    coefs <- as.matrix((first) \%*\% t(X) \%*\% y)
    resids <- (diag(nrow(X)) - (X \%*\% first \%*\% t(X))) \%*\% y
    varyest <- (t(resids) \%*\% resids)/(nrow(X) - ncol(X))
    mylist <- list(coefficients = coefs, vcov = first * as.numeric(varyest))
    return(mylist)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
