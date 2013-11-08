#define function
fastlm<-function(X, y, na.rm = FALSE) 
  {
    n<-length(y)
    p<-ncol(X)
    
    # calculating transpose(X)%*%X
    A<-crossprod(X,X) 
    C<-crossprod(X,y)
    
    #cholesky decomposition
    Q<-chol(A)
   
    #solve betahat
    temp1<-forwardsolve(t(Q),C) 
    betahat<-backsolve(Q,temp1) 
    
    #calculate covirance of beta
    cov_beta<-chol2inv(Q)*as.numeric(crossprod(y-X%*%betahat)/(n-p))
    
    return(list(coeffients=betahat,vcov=cov_beta))
  }
set.seed(2)
## Generate predictor matrix
n <- 100000
p <- 500
X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))

## Coefficents
b <- rnorm(p)

## Response
y <- X %*% b + rnorm(n)


fit <- fastlm(X, y)
str(fit)
system.time(fit<-fastlm(X, y))