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
### Title: Fast Multivariate Normal Density
### Aliases: dmvnorm mvn.density

### ** Examples

        ## Create the covariance matrix
        n <- 100
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
        
        ## The function is currently defined as
        dmvnorm <- function(x, mu, S, log = TRUE) {
        if (!is.matrix(x)){
                x=t(as.matrix(x))
        }
        k=length(mu)
        n=nrow(x)
        
        ##Check if S is positive definite
        R=tryCatch({chol(S)},
                 error=function(e){
                         message("S is not positive definite")
                 })
        
        logdetS=2*sum(log(diag(R)))
        T=x-rep(1,n)
        C=forwardsolve(t(R),t(T))
        term3=diag(crossprod(C))
        fx=(-k/2)*log(2*pi)-(1/2)*logdetS-(1/2)*term3
        if(log==TRUE){
                fx=fx
        }else {
                fx=exp(fx)
        }
        return(fx)
        }



cleanEx()
nameEx("fastlm")
### * fastlm

flush(stderr()); flush(stdout())

### Name: fastlm
### Title: Fast Linear Model Fitting
### Aliases: fastlm fast.linear

### ** Examples

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
        
        
        ## The function is currently defined as
        fastlm=function(X,y,na.rm=FALSE){
        if(na.rm==TRUE){
                r=cbind(X,y)
                X=X[complete.cases(r),]
                y=as.matrix(y[complete.cases(r)])
        }
        
        ##Cholesky factorization for coefficients
        A=crossprod(X)
        B=crossprod(X,y)
        R=chol(A)
        Rbeta=forwardsolve(t(R),B)
        coefficients=backsolve(R,Rbeta)
        
        ##Calculate VCOV
        n=length(y)
        p=ncol(X)
        sigmahat2=(crossprod(y)-crossprod(coefficients,B))/(n-p)
        Ainv=chol2inv(R)
        vcov=as.numeric(sigmahat2)*Ainv
        
        return(list(coefficients=coefficients,vcov=vcov))
        }



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
