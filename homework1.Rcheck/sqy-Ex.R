pkgname <- "sqy"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('sqy')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("dmvnorm")
### * dmvnorm

flush(stderr()); flush(stdout())

### Name: dmvnorm
### Title: Evaluate multivariate normal density
### Aliases: dmvnorm
### Keywords: ~kwd1 ~kwd2

### ** Examples
n <- 100
n2 <- n^2
xg <- seq(0, 1, length = n)
yg <- xg
g <- data.matrix(expand.grid(xg, yg))
D <- as.matrix(dist(g))
phi <- 5

S <- exp(-phi * D)
mu <- rep(0, n2)
set.seed(1)
x <- matrix(rnorm(n2), byrow = TRUE, ncol = n2)
dmvnorm(x, mu, S, log = TRUE)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, mu, S, log = TRUE) 
{
    R = tryCatch({
        chol(S)
    }, error = function(e) {
        message("S is not positive definite")
    })
    logs = 2 * sum(log(diag(R)))
    y = x - rep(1, nrow(x)) %*% t(mu)
    z = forwardsolve(t(R), t(y))
    cov = diag(crossprod(z))
    final = (-length(mu)/2) * log(2 * pi) - (1/2) * logs - (1/2) * 
        cov
    if (log == TRUE) {
        final = final
    }
    else {
        final = exp(final)
    }
    return(final)
  }



cleanEx()
nameEx("fastlm")
### * fastlm

flush(stderr()); flush(stdout())

### Name: fastlm
### Title: Conduct a fast linear model regression
### Aliases: fastlm
### Keywords: ~kwd1 ~kwd2

### ** Examples
set.seed(2)
n <- 100000
p <- 500
X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
b <- rnorm(p)
y <- X 
fit <- fastlm(X, y)
str(fit)

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, y, na.rm = FALSE) 
{
    if (na.rm == TRUE) {
        cb = cbind(x, y)
        X = X[complete.cases(cb), ]
        y = as.matrix(y[complete.cases(cb)])
    }
    n = nrow(X)
    p = ncol(X)
    l = chol(crossprod(X, X))
    cp = crossprod(X, y)
    cof = forwardsolve(t(l), forwardsolve(t(l), cp), transp = TRUE)
    vcov <- as.numeric((crossprod(y, y) - crossprod(cof, cp))) * 
        chol2inv(l)/(n - p)
    list(cof, vcov)
  }



cleanEx()
nameEx("sqy-package")
### * sqy-package

flush(stderr()); flush(stdout())

### Name: sqy-package
### Title: What the package does (short line) ~~ package title ~~
### Aliases: sqy-package sqy
### Keywords: package

### ** Examples

n <- 100
n2 <- n^2
xg <- seq(0, 1, length = n)
yg <- xg
g <- data.matrix(expand.grid(xg, yg))
D <- as.matrix(dist(g))
phi <- 5

S <- exp(-phi * D)
mu <- rep(0, n2)
set.seed(1)
x <- matrix(rnorm(n2), byrow = TRUE, ncol = n2)
dmvnorm(x, mu, S, log = TRUE)
set.seed(2)
n <- 100000
p <- 500
X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
b <- rnorm(p)
y <- X 
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
