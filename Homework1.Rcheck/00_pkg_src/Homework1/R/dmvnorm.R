dmvnorm <-
function(x, mu, S, log = TRUE) {
    ## x = data vector, mu = mean, S = covariance matrix
    ## need: fast quadratic form
    ## fast determinant of positive definite matrix

    ## if x numeric array, make row vector
    if(!is.matrix(x)) x <- matrix(x, nrow=1, ncol=ncol(S))
    k <- ncol(S)

    ## check if positive definite
    ## cheat by having chol check
    U <- try(chol(S), silent=TRUE)
    if(class(U) == "try-error") stop("S is not positive definite")

    d <- diag(U)
    ## determinant is squared product of diagonals from decomposition.
    ## Do on log scale to avoid underflow.
    logd <- sum(log(d))

    ## quickly evaluate quadratic form
    b <- crossprod(forwardsolve(t(U), t(x-mu)))
    if(is.matrix(b)) b <- diag(b)
    logf  <- -k/2 * log(2*pi) - logd - 0.5 * b
    if(log) return(logf)
    else return(exp(logf))
}
