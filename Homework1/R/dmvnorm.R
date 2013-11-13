
dmvnorm <- function(x, mu, S, log = TRUE) {
      
      
      if (is.vector(mu)==F)  mu<-as.vector(mu)
      if (is.matrix(x)==F)  x<-as.matrix(x)
      
      if (min(eigen(S)$values)>0){
            k <- length(mu)
        rooti <- backsolve(chol(S),diag(k))
        quads <- colSums((crossprod(rooti,(x-mu)))^2)
        y=exp(-(k/2)*log(2*pi) + sum(log(diag(rooti))) - .5*quads)

        if(log==TRUE) return (log(y))
        else return (y)
      }

        else
                stop ("S is not positive definite")
        
        
}


##  Reference: http://gallery.rcpp.org/articles/dmvnorm_arma/