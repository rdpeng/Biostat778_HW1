dmvnorm = function(x, mu, S, log = TRUE) {
     if (!is.matrix(x)) {
          x = t(as.matrix(x))
     }
     n = nrow(x)
     k = ncol(x)
     
     # Check that S is positive definite
#      eigenvals = eigen(S,symmetric = TRUE, only.values = TRUE)$values
#      if (sum(eigenvals < 0) > 0) # Has negative eigenvalues
#           stop("S is not positive definite")
     
     # Use Cholesky decomposition to invert S
     tryCatch( {
          cholS = chol(S)
     }, error = function(err) {
          stop("S is not positive definite")
     })
#      cholSinv = solve(cholS)
#      Sinv = tcrossprod(cholSinv)
     z = forwardsolve(t(cholS),diag(nrow(S)))
     Sinv = backsolve(cholS,z)
     
     # Calculate the 3 terms in the exponent of the density
     term1 = -k*log(2*pi)/2
     term2 = -.5*log(prod(diag(cholS))^2) # Calculate det(S) using dets of triangular matrices
     
     # Calculate density at each point
     density = rep(0,n)
     for (i in 1:n) {
          curr_row = x[i,]
          term3 = -.5*crossprod((curr_row-mu),Sinv%*%(curr_row-mu))
          density[i] = exp(term1+term2+term3)
     }
     
     if (log) {
          logd = log(density)
          return(logd)
     }
     else
          return(density)
}