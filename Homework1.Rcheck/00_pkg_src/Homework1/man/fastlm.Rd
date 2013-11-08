\name{fastlm}
\alias{fastlm}
\alias{fast.linear}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
        Fast Linear Model Fitting
%%  ~~function to do ... ~~
}
\description{
        fit a linear regression model to outcome data y and predictor data in a matrix X.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
        fastlm(X, y, na.rm = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
        X is a design matrix including predictor data and intercept.
%%     ~~Describe \code{X} here~~
}
  \item{y}{
        y is a vector of outcome data.
%%     ~~Describe \code{y} here~~
}
  \item{na.rm}{
        na.rm=TRUE to remove NA values otherwise the default value is FALSE.
%%     ~~Describe \code{na.rm} here~~
}
}
\details{
This function fits a linear regression model to outcome data y and predictor data in a matrix X. Matrix X also includes the intercept,i.e. the design matrix. This is faster than lm.fit().
%%  ~~ If necessary, more details than the description above ~~
}
\value{
        \item{coefficients}{
                The coefficients returned are the fitted regression coefficient to the input dataset.
        }
        \item{vcov}{
                The vcov returned is the variance-covariance matrix of the estimated coefficients.
        }
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Advanced Statistics Computing Class, Dr. Peng.
}
\author{
        Yu Du
%%  ~~who you are~~
}
\examples{
        set.seed(2)
        ## Generate predictor matrix
        n <- 100000
        p <- 500
        X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
        
        ## Coefficents
        b <- rnorm(p)
        
        ## Response
        y <- X \%*\% b + rnorm(n)
        
        fit <- fastlm(X, y)
        str(fit)
        
        
        ## The function is currently defined as
        function (X, y, na.rm = FALSE) 
        {
                if (na.rm == TRUE) {
                        r = cbind(X, y)
                        X = X[complete.cases(r), ]
                        y = as.matrix(y[complete.cases(r), ])
                }
                A = crossprod(X)
                B = crossprod(X, y)
                R = chol(A)
                Rbeta = forwardsolve(t(R), B)
                coefficients = backsolve(R, Rbeta)
                n = length(y)
                p = ncol(X)
                sigmahat2 = crossprod(y - X \%*\% coefficients)/(n - p)
                Ainv = chol2inv(R)
                vcov = as.numeric(sigmahat2) * Ainv
                return(list(coefficients = coefficients, vcov = vcov))
        }
        
        
}
