\name{dmvnorm}
\alias{dmvnorm}
\alias{mvn.density}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
        Fast Multivariate Normal Density
%%  ~~function to do ... ~~
}
\description{
Evaluates the k-dimensional multivariate Normal density with mean mu and covariance matrix S.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dmvnorm(x, mu, S, log = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
        x is a n-by-k matrix with each row is a sample from a k-dimensional multivariate normal distribution
%%     ~~Describe \code{x} here~~
}
  \item{mu}{
        mu is the mean vector of length k for the given multivariate normal distribution.
%%     ~~Describe \code{mu} here~~
}
  \item{S}{
        S is the p-by-p variance covariance matrix for the given multivariate normal distribution.
%%     ~~Describe \code{S} here~~
}
  \item{log}{
        if log=TRUE which is the default value, the logged density values will be returned otherwise the original density value will be returned.
%%     ~~Describe \code{log} here~~
}
}
\details{
This function evaluates the k-dimensional multivariate Normal density with mean mu and covariance matrix S.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A vector of length n of density values will be returned. n is the number of rows of input matrix x.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
        Advanced Statistics Computing Class, Dr. Peng.
%% ~put references to the literature/web site here ~
}
\author{
        Yu Du           
%%  ~~who you are~~
}
\examples{
        ## Create the covariance matrix
        n <- 100
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
        
        ## The function is currently defined as
        function (x, mu, S, log = TRUE) 
        {
                k = length(mu)
                n = nrow(x)
                R = tryCatch({
                        chol(S)
                }, error = function(e) {
                        message("S is not positive definite")
                })
                logdetS = 2 * sum(log(diag(R)))
                T = x - rep(1, n) \%*\% t(mu)
                C = forwardsolve(t(R), t(T))
                term3 = diag(crossprod(C))
                fx = (-k/2) * log(2 * pi) - (1/2) * logdetS - (1/2) * term3
                if (log == TRUE) {
                        fx = fx
                }
                else {
                        fx = exp(fx)
                }
                return(fx)
        }
}

