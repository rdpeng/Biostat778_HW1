dmvnorm<-function(x, mu, S, log=TRUE)
{
      cs<-tryCatch({chol(S)}, error=function(cond)
            {stop("S is not positive definite")})
      
      ## calculate log density
      a<-forwardsolve(t(cs),t(x)-mu)
      logf<--log(2*pi)*length(mu)/2-sum(log(diag(cs)))-colSums(a*a)/2

      if(log==T) result<-logf
      else result<-exp(logf)
      result
}