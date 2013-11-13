\name{fastlm}
\alias{fastlm}
\title{Fast Linear Regression}
\description{
Faster linear regression than basic lm() function 
}
\usage{
fastlm(X,y,na.rm=FALSE)
}
\arguments{
\item{X}{The covariate matrix}
\item{y}{The response variable}
\item{na.rm}{Whether observations with na should be excluded}
}
\details{
Make sure that X is a covariate matrix
}
\value{
A list of coefficients and covariance matrix
}