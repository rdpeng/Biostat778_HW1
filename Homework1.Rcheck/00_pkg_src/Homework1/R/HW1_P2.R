#define function dmvnorm
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
               message("S is not positive definite")
             })
  
  #compute Q_inverse*(x-mu) 
  #note that t(x-mu)%*%inv(Q%*%t(Q))%*%(x-mu)=crossprod(inv(Q)%*%(x-mu))
  #the easiest way to compute that is to solve t(Q)%*%(inv(Q)%*%(x-mu))=x-mu
        A=forwardsolve(t(Q),t(x)-mu)
  
  #compute density
        density=(-k/2)*log(2*pi)-sum(log(diag(Q)))-(1/2)*diag(crossprod(A))
  
  #check if log argument
        if(log!=TRUE){
            density=exp(density)
        }  
        return(density)
}