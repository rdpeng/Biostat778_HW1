##Fast Multivariate Normal Density
dmvnorm <- function(x, mu, S, log = TRUE) {
        k=length(mu)
        n=nrow(x)
        
        ##Check if S is positive definite
        R=tryCatch({chol(S)},
                 error=function(e){
                         message("S is not positive definite")
                 })
        
        logdetS=2*sum(log(diag(R)))
        T=x-rep(1,n)%*%t(mu)
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

