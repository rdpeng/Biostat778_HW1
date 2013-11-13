\name{fastlm}
\alias{fastlm}
\title{
fast linear regression estimation
}
\description{
compute parameter MLE and its var-cov matrix for a linear regression model
}
\usage{
fastlm(X, y, na.rm = FALSE)
}
\arguments{
  \item{X}{
an n by p input matrix
}
  \item{y}{
an n by 1 input vector 
}
  \item{na.rm}{
missing values in X or y should be removed
}
}
\value{
reture a list contains parameter MLE and its var-cov matrix for a linear regression model
}