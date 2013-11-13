dmvnorm = function(x, mu, S, log = TRUE) {
     if (!is.matrix(x)) {
          x = t(as.matrix(x))
     }
     n = nrow(x)
     k = ncol(x)
     
     # Use Cholesky decomposition to invert S
     tryCatch( {
          cholS = chol(S)
     }, error = function(err) {
          stop("S is not positive definite")
     })
#      cholSinv = solve(cholS)
#      Sinv = tcrossprod(cholSinv)
#      z = forwardsolve(t(cholS),diag(nrow(S)))
#      Sinv = backsolve(cholS,z)
     
     # Calculate the 3 terms in the exponent of the density
     term1 = -k*log(2*pi)/2
     term2 = -.5*log(prod(diag(cholS))^2) # Calculate det(S) using dets of triangular matrices
     
     # Split exponent into exp(-.5QQ')
     # Q = (X-mu)'U^-1 where U = chol(S)
#      Q = forwardsolve(t(cholS),t(x)-mu)
     Q = forwardsolve(cholS,t(x)-mu,transpose=TRUE,upper.tri=TRUE)
     exp_term = -.5*colSums(Q*Q)
#      exp_term = -.5*crossprod(Q)
     
     # Calculate normalizing constants
#      coef1 = (2*pi)^(-k/2)
#      coef2 = 1/sqrt(prod(diag(cholS)))
     
     # Calculate density
#      density = coef1*coef2*exp(exp_term)
     density = exp(term1+term2+exp_term)
     
     # Calculate density at each point
#      density = rep(0,n)
#      for (i in 1:n) {
#           curr_row = x[i,]
#           term3 = -.5*crossprod((curr_row-mu),Sinv%*%(curr_row-mu))
#           density[i] = exp(term1+term2+term3)
#      }
     
     if (log)
          return(log(density))
     else
          return(density)
}