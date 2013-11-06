dmvnorm <-
function(x, mu, S, log = TRUE) {
  U = chol(S)                              # by using cholesky decomposition, the positive-definity is checked.
  log_detS_sqrt = sum(log(diag(U)))        # compute the 1/2 log (det(S)), avoid computing product of many small numbers 
  k = dim(x)[2]
  fi = function(xi){                       # define the function for evaluating the density
    D=xi-mu
    D1 = forwardsolve(t(U),D)              # compute the s^{-1/2}*(x-mu), which is a colum vector
    D2 = sum(D1^2)                         # compute the mahalanobis distance
    -k/2*log(2*pi)-log_detS_sqrt-1/2*D2    # the log density
  }
  if (log==TRUE){
    mvpdf = apply(x,1,fi)                  # apply the function to each row of x
  }
  else{                                    # the unlogged density
    mvpdf = exp(apply(x,1,fi))
  }
  mvpdf                                    # return density
}
