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
### Title: Evaluate the density of a multivariate normal distribution
### Aliases: dmvnorm

### ** Examples

x <- matrix(rnorm(10*9), ncol=9)
mu <- rep(0,9)
xg <- seq(0, 1, length = 3)
yg <- xg
g <- data.matrix(expand.grid(xg, yg))
D <- as.matrix(dist(g))
S <- exp(D * -1)
dmvnorm(x, mu, S)



cleanEx()
nameEx("fastlm")
### * fastlm

flush(stderr()); flush(stdout())

### Name: fastlm
### Title: Fast linear regression
### Aliases: fastlm

### ** Examples

X <- matrix(rnorm(10),ncol=2)
y <- rnorm(5)
fastlm(X, y)



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
