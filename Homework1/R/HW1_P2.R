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
        A=forwardsolve(t(Q),t(temp1))
        temp2=diag(crossprod(A))
        
        #Computing density
        density=(-k/2)*log(2*pi)-(1/2)*2*sum(log(diag(Q)))-(1/2)*temp2
        
        #Check if log is false
        if(log==FALSE){
                density=exp(density)
        }
        
        return(density)
}

