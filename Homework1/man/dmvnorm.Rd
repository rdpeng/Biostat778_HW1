\name{dmvnorm}
\alias{dmvnorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Computing the Multivariate Normal Density MORE Efficiently
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function dmvnorm() evaluates the k-dimensional multivariate Normal density with mean mu and covariance S
}
\usage{
dmvnorm(x, mu, S, log = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
a n*k matrix of points to be evaluated}
  \item{mu}{
%%     ~~Describe \code{mu} here~~
a length k vector of means for the k-dimensional Normal}
  \item{S}{
%%     ~~Describe \code{S} here~~
a k*k covariance matrix}
  \item{log}{
%%     ~~Describe \code{log} here~~
whether the return should be the log density. Default is TRUE}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
This funtion returns a vector of length n containing the values of the multivariate Normal density evaluated at the n points. If log = TRUE, it returns the log density at those points. If log is NOT TRUE, it will return the density values.}
\references{
%% ~put references to the literature/web site here ~
PH140.778 Advanced Statistical Computing Dr.Peng
}
\author{
%%  ~~who you are~~
Lu Li}

\examples{
%The function is currently defined as
dmvnorm <- function(x, mu, S, log=TRUE) {
  #k dimension multivariate normal
        k=length(mu) 
        if(is.matrix(x)==FALSE){
        x=as.matrix(t(x))
        }
  #n data points
        n=nrow(x)
  
  #check positive definite by trying to Cholesky decomposition
        Q=tryCatch({chol(S)},
             error=function(li){
               message("S cannot be a covariance matrix")
             })
  
  #compute Q_inverse*(x-mu) 
  #note that t(x-mu)\%*\%inv(Q\%*\%t(Q))\%*\%(x-mu)=crossprod(inv(Q)\%*\%(x-mu))
  #the easiest way to compute that is to solve t(Q)\%*\%(inv(Q)\%*\%(x-mu))=x-mu
        temp1=x-rep(1,n)\%*\%t(mu)
        A=forwardsolve(t(Q),t(temp1))
        temp2=diag(crossprod(A))
  
  #compute density
        density=(-k/2)*log(2*pi)-(1/2)*2*sum(log(diag(Q)))-(1/2)*temp2
  
  #check if log argument
        if(log!=TRUE){
            density=exp(density)
        }  
        return(density)
}
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

dmvnorm(x, mu, S, log = TRUE)

}