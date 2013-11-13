dmvnorm <-
function(x, mu, S, log = TRUE) {
  ##Check if S is positive definite
  R=tryCatch({chol(S)},
             error=function(e){
               message("S is not positive definite")
             })
  
  logs=2*sum(log(diag(R)))
  y=x-rep(1,nrow(x))%*%t(mu)
  z=forwardsolve(t(R),t(y))
  cov=diag(crossprod(z))
  final=(-length(mu)/2)*log(2*pi)-(1/2)*logs-(1/2)*cov
  if(log==TRUE){
    final=final
  }else {
    final=exp(final)
  }
  return(final)
}
