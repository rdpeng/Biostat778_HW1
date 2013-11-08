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
### Title: Fast Multivariate Normal (log)density Calculation
### Aliases: dmvnorm
### Keywords: multivariate normal,density,fast

### ** Examples

## Create the covariance matrix
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



cleanEx()
nameEx("fastlm")
### * fastlm

flush(stderr()); flush(stdout())

### Name: fastlm
### Title: Fast Linear Model Fit
### Aliases: fastlm
### Keywords: linear model,multivariate normal,fast

### ** Examples

set.seed(2)
## Generate predictor matrix
n <- 100000
p <- 500
X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))

## Coefficents
b <- rnorm(p)

## Response
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
