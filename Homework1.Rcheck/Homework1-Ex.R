pkgname <- "Homework1"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('Homework1')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("dmvnorm")
### * dmvnorm

flush(stderr()); flush(stdout())

### Name: dmvnorm
### Title: multivariate normal density
### Aliases: dmvnorm
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, mu, S, log = TRUE) 
{
    if (!is.matrix(x)) 
        x <- matrix(x, nrow = 1, ncol = ncol(S))
    k <- ncol(S)
    U <- try(chol(S), silent = TRUE)
    if (class(U) == "try-error") 
        stop("S is not positive definite")
    d <- diag(U)
    logd <- sum(log(d))
    b <- crossprod(forwardsolve(t(U), t(x - mu)))
    if (is.matrix(b)) 
        b <- diag(b)
    logf <- -k/2 * log(2 * pi) - logd - 0.5 * b
    if (log) 
        return(logf)
    else return(exp(logf))
  }



cleanEx()
nameEx("fastlm")
### * fastlm

flush(stderr()); flush(stdout())

### Name: fastlm
### Title: Fast linear model
### Aliases: fastlm
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, y, na.rm = FALSE) 
{
    if (na.rm) {
        if (any(is.na(X))) {
            inx <- which(is.na(rowSums(X)))
            X <- X[-inx, , drop = FALSE]
            y <- y[-inx]
        }
        if (any(is.na(y))) {
            iny <- which(is.na(y))
            y <- y[-iny]
            X <- X[-iny, , drop = FALSE]
        }
    }
    n <- nrow(X)
    p <- ncol(X)
    U <- chol(t(X) %*% X)
    b <- backsolve(U, forwardsolve(t(U), t(X) %*% y))
    sigma2 <- 1/(n - p) * (t(y) %*% y - t(b) %*% t(X) %*% y)
    s2 <- diag(as.numeric(sigma2), p)
    var.b <- backsolve(U, forwardsolve(t(U), s2))
    return(list(coefficients = b, vcov = var.b))
  }



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
