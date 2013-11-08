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


