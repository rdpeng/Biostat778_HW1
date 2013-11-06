fastlm <-
function(X,y,na.rm=FALSE){ # using cholesky decomposition.
  size = dim(X)
  n = size[1]
  p = size[2]
  A = t(X)%*%X                    # compute X'X ; complexity is about np^2
  B = t(X)%*%y                    # compute X'Y
  U = chol(A)                     # apply cholesky decomposition : X'X = U'U, U is upper triangle; complexity is about (p^3)/3
  y1 = forwardsolve(t(U),B)       # solve U'y = X'Y, y = U beta
  beta = backsolve(U,y1)          # solve U beta = y 
  XtXinv = chol2inv(U)            # compute (X'X)^{-1}
  resid = y-X%*%beta              # compute the residuals
  sigma2 = sum(resid^2)/(n-p)     # estimate the variance of Y: Y ~ N(X * beta, sigma2 * I)
  vcov = sigma2*XtXinv            # Var(beta_hat) = sigma2 * (X'X)^{-1}
  list(coefficients=beta, vcov=vcov)
}
