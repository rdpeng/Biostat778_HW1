##########################################
## A faster linear regression.
## Arguments:
##      X = an nxp matrix
##      y = a vector of length n
##      na.rm
##########################################
fastlm <- function(X, y, na.rm=FALSE) {
    ## take advantage of X'X being symmetric
    if(na.rm) {
        if(any(is.na(X))) {
            inx <- which(is.na(rowSums(X)))
            X <- X[-inx,, drop=FALSE]
            y <- y[-inx]
        }
        if(any(is.na(y))) {
            iny <- which(is.na(y))
            y <- y[-iny]
            X <- X[-iny,, drop=FALSE]
        }
    }

    n <- nrow(X)
    p <- ncol(X)
    U <- chol(t(X)%*%X)

    b <- backsolve(U, forwardsolve(t(U), t(X)%*%y))
    sigma2 <- 1/(n-p) * (t(y) %*% y - t(b) %*% t(X) %*% y)
    s2 <- diag(as.numeric(sigma2), p)
    var.b <- backsolve(U, forwardsolve(t(U), s2))
    return(list("coefficients" = b, "vcov" = var.b))
}
