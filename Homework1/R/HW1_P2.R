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
               message("S cannot be a covariance matrix")
             })
  
  #compute Q_inverse*(x-mu) 
  #note that t(x-mu)%*%inv(Q%*%t(Q))%*%(x-mu)=crossprod(inv(Q)%*%(x-mu))
  #the easiest way to compute that is to solve t(Q)%*%(inv(Q)%*%(x-mu))=x-mu
        temp1=x-rep(1,n)%*%t(mu)
        A=forwardsolve(t(Q),t(temp1))
        temp2=diag(crossprod(A))
  
  #compute density
        density=(-k/2)*log(2*pi)-(1/2)*2*sum(log(diag(Q)))-(1/2)*temp2
  
  #check if log argument
        if(log!=TRUE){
            density=exp(density)
        }  
        return(density)
}