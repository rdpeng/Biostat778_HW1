pkgname <- "Homework1"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('Homework1')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("dmvnorm")
### * dmvnorm

flush(stderr()); flush(stdout())

### Name: dmvnorm
### Title: Computing the Multivariate Normal Density MORE Efficiently
### Aliases: dmvnorm

### ** Examples



## The function is currently defined as
function (x, mu, S, log = TRUE) 
{
    k = length(mu)
    n = nrow(x)
    Q = tryCatch({
        chol(S)
    }, error = function(li) {
        message("S cannot be a covariance matrix")
    })
    temp1 = x - rep(1, n) %*% t(mu)
    A = forwardsolve(t(Q), t(temp1))
    temp2 = diag(crossprod(A))
    density = (-k/2) * log(2 * pi) - (1/2) * 2 * sum(log(diag(Q))) - 
        (1/2) * temp2
    if (log == FALSE) {
        density = exp(density)
    }
    return(density)
  }



cleanEx()
nameEx("fastlm")
### * fastlm

flush(stderr()); flush(stdout())

### Name: fastlm
### Title: Faster Way to Fit Linear Regression Models
### Aliases: fastlm
### Keywords: Cholesky decomposition Linear Regression Model

### ** Examples

function (X, y, na.rm = FALSE) 
{
    n <- length(y)
    p <- ncol(X)
    
    ##Check if missing values in X and y should be removed
    if (na.rm == TRUE) {
        Z = cbind(X, y)
        X = X[complete.cases(Z), ]
        y = as.matrix(y[complete.cases(Z)])
    }
    A <- crossprod(X)
    C <- crossprod(X, y)
    
    ##Cholesky decomposition
    Q <- chol(A)
    temp1 <- forwardsolve(t(Q), C)
    betahat <- backsolve(Q, temp1)
    cov_beta <- chol2inv(Q) * as.numeric(crossprod(y - X %*% 
        betahat)/(n - p))
    return(list(coeffients = betahat, vcov = cov_beta))
  }
    set.seed(2)
## Generate predictor matrix
    n <- 100
    p <- 5
    X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))

## Coefficents
    b <- rnorm(p)

## Response
    y <-X%*%b + rnorm(n)

    fit <- fastlm(X, y)
    str(fit)



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
