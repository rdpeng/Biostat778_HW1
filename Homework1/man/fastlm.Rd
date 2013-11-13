\name{fastlm}
\alias{fastlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Conduct a fast linear model regression
%%  ~~function to do ... ~~
}
\description{The function can conduct a fast linear model regression.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fastlm(X, y, na.rm = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a n×p matrix of variables
%%     ~~Describe \code{X} here~~
}
  \item{y}{a vector of observations
%%     ~~Describe \code{y} here~~
}
  \item{na.rm}{deciding whether NA should be removed
%%     ~~Describe \code{na.rm} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{The function returns a list of 2:coefficients and vcov
	
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{JerryPPA
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{set.seed(2)
n <- 100000
p <- 500
X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
b <- rnorm(p)
y <- X %*% b + rnorm(n)
source("fastlm.R")
fit <- fastlm(X, y)
str(fit)

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, y, na.rm = FALSE) 
{
    if (na.rm == TRUE) {
        cb = cbind(x, y)
        X = X[complete.cases(cb), ]
        y = as.matrix(y[complete.cases(cb)])
    }
    n = nrow(X)
    p = ncol(X)
    l = chol(crossprod(X, X))
    cp = crossprod(X, y)
    cof = forwardsolve(t(l), forwardsolve(t(l), cp), transp = TRUE)
    vcov <- as.numeric((crossprod(y, y) - crossprod(cof, cp))) * 
        chol2inv(l)/(n - p)
    list(cof, vcov)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
