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

