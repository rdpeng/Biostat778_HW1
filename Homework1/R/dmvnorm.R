###Fast Multivariate Normal Density
dmvnorm <- function(x, mu, S, log=TRUE) {
      if (is.vector(x))
            x <- as.matrix(t(x))
      cholS <- tryCatch(chol(S),error=function(cond) {stop("S is not positive definite",call. = F)})
      k <- ncol(x)
      res <- sapply(1:nrow(x), function(i) {
            -(k/2)*log(2*pi) - sum(log(diag(cholS))) - .5*sum(crossprod(forwardsolve(t(cholS),x[i,]-mu)))
      })
      if (log == TRUE) {
            res
      } else {
            exp(res)
      }
}

