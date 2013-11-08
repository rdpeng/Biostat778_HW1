#Define function Fast Multivariate Normal Density
dmvnorm <- function(x, mu, S, log = TRUE) {
        k=length(mu) 
        n=nrow(x)
        
        #Check positive definite
        Q=tryCatch({chol(S)},
                   error=function(li){
                           message("S cannot be a covariance matrix")
                   })
        
        #Computing t(x-mu)*Q_inverse
        temp1=x-rep(1,n)%*%t(mu)
        C=forwardsolve(t(Q),t(temp1))
        temp2=diag(crossprod(C))
        logdensity=(-k/2)*log(2*pi)-(1/2)*2*sum(log(diag(Q)))-(1/2)*temp2
        if(log==FALSE){
                logdensity=exp(logdensity)
        }
        
        return(logdensity)
}

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