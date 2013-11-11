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
### Title: Fast multivariate normal density
### Aliases: dmvnorm

### ** Examples

## Create the covariance matrix
n <- 10
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



cleanEx()
nameEx("fastlm")
### * fastlm

flush(stderr()); flush(stdout())

### Name: fastlm
### Title: Fast linear regression
### Aliases: fastlm

### ** Examples

set.seed(2)
## Generate predictor matrix
n <- 1000
p <- 50
X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))

## Coefficents
b <- rnorm(p)

## Response
y <- X 

fit <- fastlm(X, y)



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
