dmvnorm<-function(x, mu, S, log=TRUE)
{
      k<-length(mu)
      cs<-tryCatch({chol(S)}, error=function(cond)
            {stop("S is not positive definite")})
      
      ## calculate log density
      a<-forwardsolve(t(cs),t(x)-mu)
      logf<--log(2*pi)*k/2-sum(log(diag(cs)))-colSums(a^2)/2

      if(log) result<-logf
      else result<-exp(logf)
      result
}