###Fast Linear Regression
fastlm <- function(X, y, na.rm = FALSE) {
      if (na.rm == T) {
            nonnaid <- complete.cases(X,y)
            X <- X[nonnaid,]
            y <- y[nonnaid]
      }
      k <- ncol(X)
      XtXinv <- chol2inv(chol(crossprod(X)))
      Xty <- crossprod(X, y)
      coef <- XtXinv %*% Xty
      var <- (sum(y^2)-sum(Xty * coef))/ (length(y)-k) * XtXinv      
      list(coefficients=coef,vcov=var)
}
