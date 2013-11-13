#' Do fast linear models
#'
#' \code{fastlm} Does lm, but faster. Uses cholesky decomposition of crossprod(X) to invert X'X. 
#'
#' @param X a nxp matrix, where k is the number of measurements per observation, and n is the number of observations
#' @param y a n-length vector of outcomes
#' @param na.rm Whether the X or y data has na's. If TRUE, take the time to remove them. Otherwise skip it if you don't expect an error.
#' @return A vector of densities, length-n
#' @export
#'
#' @examples
#' set.seed(2)
#'  #####Generate predictor matrix
#' n <- 1000
#' p <- 5
#' X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
#' 
#' ## Coefficents
#' b <- rnorm(p)
#' 
#' ## Response
#' y <- X %*% b + rnorm(n)
#' 
#' 
#' 
#' 
#' system.time( fit_fast <- fastlm(X, y) )
#' system.time( fit_compare <- lm.fit(X, y) )
#' str(fit_fast)

fastlm <- function(X, y, na.rm = FALSE) {
  if(na.rm){
    naInd<-apply(cbind(y,X),1,function(x) any(is.na(x)))
    X<-X[!naInd,]
    y<-y[!naInd]
  }
  p<-dim(X)[2]
  n<-dim(X)[1]
  XtX<-crossprod(X) #faster than t(X)%*%X
  R<-chol(XtX) #works for any pos def symmetric matrix. In this case, it outputs the R from qr(X)=QR, as X'X=R'Q'QR=R'R
  Rinv<-backsolve(R,diag (rep(1,p))) #solving an upper triangular matrix
  XtXinv<-tcrossprod(Rinv) #transpose of crossprod
  Xty<-crossprod(X,y)
  beta<- XtXinv %*% Xty
  resid<- y- X%*%beta
  sd_e<-sum(resid^2)/(n-p)

  return(list(coefficients=beta, vcov=sd_e*XtXinv))
}




#' get multivariate normal density for a matrix of vectors.
#'
#' \code{dmvnorm} gets the multivariate density.
#'
#' @param x nxk matrix, where k is the number of measurements per observation, and n is the number of observations
#' @param mu k-length vector of means
#' @param S a kxk covariance matrix
#' @param log whether to output log denisites or not.
#' @return A vector of densities, length-n
#' @export
#' @examples
#' kSqrt <- 5 #our final density will be for a random vector of length 15
#' k <- kSqrt^2
#' xg <- seq(0, 1, length = kSqrt)
#' yg <- xg
#' g <- data.matrix(expand.grid(xg, yg))
#' D <- as.matrix(dist(g))
#' phi <- 5
#' S <- exp(-phi * D) #S ends up being k by k
#' n<-4 #four observations
#'
#' mu <- rep(0, k)
#' set.seed(1)
#' x <- matrix(rnorm(k*n), byrow = TRUE, ncol = k) #n x k, this is the usual input to dmvnorm i mvtnorm
#' 
#' system.time( fastOutput <-dmvnorm(x=x,mu=mu,S=S,log=TRUE) )
#' str(fastOutput)
dmvnorm <- function(x, mu, S, log = TRUE) {
  ## Your code here
  if(class(x)=='numeric') x<-matrix(x,nrow=1)
  n<-dim(x)[1]
  k<-dim(x)[2]

  #I lifted this from dmvnorm in the mvtnorm package, det(S) was failing. This also gives us a way to check positive definitelness. See equivalent conditions of pos def on the wikipedia page
  eigenValS <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  if(any(eigenValS<=0)) stop("S is not positive definite")

  logDetS <- sum(log(eigenValS)) #equal to the way of evaluating this below, but less stable.
  #logdet <- log(prod(eigenValS))# equal to the above, but less stable. It's more stable to take sums of logs than logs of products.


  #system.time({
  muMat<-matrix(rep(mu,times=n),nrow=n,ncol=k)
  xTilde<-x-muMat #n x k
  #we want S^inv X = Q, Q is the thing we want.
  #X=SQ; X and S are known. Here, X is the centered x that this func takes as an argument.
  SigmaInvXTilde<-solve(a=S,b=t(xTilde)) #k x n
  quadTerms <-rep(NA,n)
  for (i in 1:n) quadTerms[i]<-sum((SigmaInvXTilde[,i]*xTilde[i,]))
  #})#end system.time
  #quadTerms faster than dmvnorms malananobis() way to calculate this
  #system.time(distval <- mahalanobis(x, center = mu, cov = S))

  logdnorm <- -(k/2)*log(2*pi)-.5*logDetS-.5*quadTerms
  
  if(log) return(logdnorm)
  if(!log) return(exp(logdnorm))
   
}




#This alternate method uses cholesky(S) to evaluate the quadratic
#It's not as fast as the first attempt but still faster than mvtnorm package
dmvnorm_cholesky <- function(x, mu, S, log = TRUE) {
  ## Your code here
  n<-dim(x)[1]
  k<-dim(x)[2]

  #I lifted this from dmvnorm in the mvtnorm package, det(S) was failing. This also gives us a way to check positive definitelness. See equivalent conditions of pos def on the wikipedia page
  eigenValS <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  if(any(eigenValS<=0)) stop("S is not positive definite")

  logDetS <- sum(log(eigenValS)) #equal to the way of evaluating this below, but less stable.
  #logdet <- log(prod(eigenValS))# equal to the above, but less stable. It's more stable to take sums of logs than logs of products.


  #system.time({
  muMat<-matrix(rep(mu,times=n),nrow=n,ncol=k)
  xTilde<-x-muMat #n x k
  R<- chol(S) #let S:=R'R
  Rinv<- backsolve(R,diag (rep(1,k))) #for upper triangulars
  #Sinv=(Rinv)(Rinv)'
  #(x-mu)'Sinv(x-mu)
  #=crossprod((Rinv)'(x-mu))
  #=crossprod((Rinv)'(xTilde)')
  #=crossprod((xTilde%*%Rinv)')
  #=tcrossprod((xTilde%*%Rinv))
  quadTerms <-rep(NA,n)
  for (i in 1:n) quadTerms[i]<-tcrossprod((xTilde[i,]%*%Rinv))

  logdnorm <- -(k/2)*log(2*pi)-.5*logDetS-.5*quadTerms
  
  if(log) return(logdnorm)
  if(!log) return(exp(logdnorm))
}

