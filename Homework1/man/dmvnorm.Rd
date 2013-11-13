\name{dmvnorm}
\alias{dmvnorm}
\title{
calculate multivariate normal densities
}
\description{
  dmvnorm is used to fast calculate the multivariate normal densities
}
\usage{
dmvnorm(x, mu, S, log = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  	a n*k matrix of points to be evaluated. 
}
  \item{mu}{
  	a length k vector of means for the k-dimensional Normal
}
  \item{S}{
  	a k*k covariance matrix
}
  \item{log}{
    A logical variable. If log =TRUE then it returns the log density at those points. If log = FALSE then it returns the density values.
}
}
\details{
This algrithm uses the Choleski Decmoposition of matrix, with the function forwardsolve(), it creates a faster way to calculate the multivariate normal densities.
}
\value{
	Return a vector of length n containing the values of the multivariate Normal density evaluated at the n points.
}

\author{
	Yuan He
}
\note{
S should be a positive definite matrix.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
forwardsolve(), chol()
}
\examples{
	set.seed(1)
	x = matrix(rnorm(9), byrow = TRUE, ncol = 9)
	mu = rep(0, 9)
	xg <- seq(0, 1, length = 3)
	yg <- xg
	g <- data.matrix(expand.grid(xg, yg))
	D <- as.matrix(dist(g))
	phi <- 5
	S <- exp(-phi * D)
	dmvnorm(x, mu, S, log = TRUE)
## The function is currently defined as
dmvnorm <- function(x, mu, S, log=TRUE)
{
  # check if S is positive definite
  l = tryCatch(chol(S),error = function(e){message ("S is not positive definite")})
  
  # log(det(S))
  logdS = 2*sum(log(diag(l)))
  
  # (x-mu)
  n = nrow(x)
  Mx = x - rep(1,n)%*%t(mu)
  
  # Y = backsolve(l,Mx)
  Y = forwardsolve(t(l),t(Mx))
  Mx3 = diag(crossprod(Y))
  
  k = length(mu)
  fx = -(k/2)*log(2*pi)-(1/2)*logdS-(1/2)*Mx3
  
  if (log == FALSE)
    fx = exp(fx)
  else fx= fx
  
  return(fx)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariant}
\keyword{normal distribution}% __ONLY ONE__ keyword per line
