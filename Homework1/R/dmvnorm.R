###Fast Multivariate Normal Density
dmvnorm <- function(x, mu, S, log=TRUE) {
      if (is.vector(x))
            x <- as.matrix(t(x))
      cholS <- tryCatch(chol(S),error=function(cond) {stop("S is not positive definite",call. = F)})
      k <- ncol(x)
      invS <- backsolve(cholS,diag(k))
      res <- sapply(1:nrow(x), function(i) {
            -(k/2)*log(2*pi) + sum(log(diag(invS))) - .5*sum((crossprod(invS,x[i,]-mu))^2)
      })
      if (log == TRUE) {
            res
      } else {
            exp(res)
      }
}
