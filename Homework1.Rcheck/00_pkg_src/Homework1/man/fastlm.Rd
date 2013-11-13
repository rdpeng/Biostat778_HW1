\name{fastlm}
\alias{fastlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Faster Way to Fit Linear Regression Models
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Comparing to the lm.fit() function in R, this fastlm() function presents a much faster way to fit a linear regression model. 
}
\usage{
fastlm(X, y, na.rm = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
predictor variable, a n*k matrix
}
  \item{y}{
%%     ~~Describe \code{y} here~~
response variable, a vector of length n
}
  \item{na.rm}{
%%     ~~Describe \code{na.rm} here~~
argument indicating that whether missing values in X or y should be removed. Default is FALSE 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
fastlm() takes the advantages of Cholesky decomposition to compute the inverse matrices, which makes a huge improvement in the efficiency of fitting linear regression models. 
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
a list of components 
\item{coefficients}{
a vector of the regression coefficients estimated using maximum likelihood
}
\item{vcov}{
the p*p covariance matrix of the estimated regression coefficients
}
}
\references{
%% ~put references to the literature/web site here ~
PH140.778 Advanced Statistical Computing Dr.Peng
}
\author{
%%  ~~who you are~~
Lu Li
}

\examples{
%## The function is currently defined as

fastlm<-function(X, y, na.rm=FALSE) {        
    #check argument na.rm
        if (na.rm!=FALSE) {
                Z=cbind(X,y)
                X=X[complete.cases(Z),]
                y=as.matrix(y)[complete.cases(Z)]
        }
        n<-length(y)
        p<-ncol(X)
    #calculating transpose(X)\%*\%X
        A<-crossprod(X)
    #calculating transpose(X)\%*\%y    
        C<-crossprod(X,y)
    
    #cholesky decomposition
        Q<-chol(A)
   
    #solve betahat
        temp1<-forwardsolve(t(Q),C) 
        betahat<-backsolve(Q,temp1) 
    
    #calculate covirance of beta
    #note that t(e)\%*\%e=t(e)\%*\%y=t(y)\%*\%y-t(y)\%*\%X\%*\%betahat
  
        cov_beta<-chol2inv(Q)*as.numeric(crossprod(y)-crossprod(X\%*\%betahat))/(n-p)
    
        return(list(coefficients=betahat,vcov=cov_beta))
}
    set.seed(2)
## Generate predictor matrix
    n <- 100
    p <- 5
    X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))

## Coefficents
    b <- rnorm(p)

## Response
    y <-X\%*\%b + rnorm(n)

    fit <- fastlm(X, y)
    str(fit)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Cholesky decomposition}
\keyword{Linear Regression Model}% __ONLY ONE__ keyword per line
