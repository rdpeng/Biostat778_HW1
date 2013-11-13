fastlm = function(X, y, na.rm = FALSE) {
     # Remove missing values
     if (na.rm) {
          good = complete.cases(X,y)
          X = as.matrix(X[good,])
          y = as.matrix(y[good])
     }
     n = nrow(X)
     p = ncol(X)
     
     # Calculate X'X and X'Y
     xtx = crossprod(X)
     xty = crossprod(X,y)
     
     # Use Cholesky decomposition to quickly solve normal equations
     chol_xtx = chol(xtx) # Upper triangular matrix
     z = forwardsolve(t(chol_xtx),xty)
     beta = backsolve(chol_xtx,z)
     
     # Use Cholesky decomposition to find (X'X)^-1
     chol_xtx_inv = solve(chol_xtx)
     xtxinv = tcrossprod(chol_xtx_inv)
     residual = y-(X%*%beta)
     var_hat = sum(residual^2)/(n-p) # scalar
     var_beta = var_hat*xtxinv
     
     beta = unname(beta)
     var_beta = unname(var_beta)
     return(list(coefficients=beta,vcov=var_beta))
}