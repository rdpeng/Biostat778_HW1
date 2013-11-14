\name{fastlm}
\alias{fastlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Fast linear model
}
\description{
    ~fast~ linear regression
}
\usage{
fastlm(X, y, na.rm = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{na.rm}{
%%     ~~Describe \code{na.rm} here~~
}
}
\details{
    Uses Cholesky transormation.
%%  ~~ If necessary, more details than the description above ~~
}
\author{
    Stephen Cristiano
}
\note{
    works pretty well.
}


\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, y, na.rm = FALSE) 
{
    if (na.rm) {
        if (any(is.na(X))) {
            inx <- which(is.na(rowSums(X)))
            X <- X[-inx, , drop = FALSE]
            y <- y[-inx]
        }
        if (any(is.na(y))) {
            iny <- which(is.na(y))
            y <- y[-iny]
            X <- X[-iny, , drop = FALSE]
        }
    }
    n <- nrow(X)
    p <- ncol(X)
    U <- chol(t(X) \%*\% X)
    b <- backsolve(U, forwardsolve(t(U), t(X) \%*\% y))
    sigma2 <- 1/(n - p) * (t(y) \%*\% y - t(b) \%*\% t(X) \%*\% y)
    s2 <- diag(as.numeric(sigma2), p)
    var.b <- backsolve(U, forwardsolve(t(U), s2))
    return(list(coefficients = b, vcov = var.b))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
