pkgname <- "Homework1"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('Homework1')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("dmvnorm")
### * dmvnorm

flush(stderr()); flush(stdout())

### Name: dmvnorm
### Title: Computing the Multivariate Normal Density MORE Efficiently
### Aliases: dmvnorm

### ** Examples

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
  #note that t(x-mu)%*%inv(Q%*%t(Q))%*%(x-mu)=crossprod(inv(Q)%*%(x-mu))
  #the easiest way to compute that is to solve t(Q)%*%(inv(Q)%*%(x-mu))=x-mu
        A=forwardsolve(t(Q),t(x)-mu)
        temp2=diag(crossprod(A))
  
  #compute density
        density=(-k/2)*log(2*pi)-sum(log(diag(Q)))-(1/2)*temp2
  
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




cleanEx()
nameEx("fastlm")
### * fastlm

flush(stderr()); flush(stdout())

### Name: fastlm
### Title: Faster Way to Fit Linear Regression Models
### Aliases: fastlm
### Keywords: Cholesky decomposition Linear Regression Model

### ** Examples


fastlm<-function(X, y, na.rm=FALSE) {        
    #check argument na.rm
        if (na.rm!=FALSE) {
                Z=cbind(X,y)
                X=X[complete.cases(Z),]
                y=as.matrix(y)[complete.cases(Z)]
        }
        n<-length(y)
        p<-ncol(X)
    #calculating transpose(X)%*%X
        A<-crossprod(X)
    #calculating transpose(X)%*%y    
        C<-crossprod(X,y)
    
    #cholesky decomposition
        Q<-chol(A)
   
    #solve betahat
        temp1<-forwardsolve(t(Q),C) 
        betahat<-backsolve(Q,temp1) 
    
    #calculate covirance of beta
    #note that t(e)%*%e=t(e)%*%y=t(y)%*%y-t(y)%*%X%*%betahat
  
        cov_beta<-chol2inv(Q)*as.numeric(crossprod(y)-crossprod(X%*%betahat))/(n-p)
    
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
    y <-X%*%b + rnorm(n)

    fit <- fastlm(X, y)
    str(fit)



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
