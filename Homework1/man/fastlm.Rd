\name{fastlm}
\alias{fastlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
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
%%  ~~ If necessary, more details than the description above ~~
}
\value{
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
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

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
