## `S' rank deficient

library(Homework1)
op <- options(scipen = 5)

set.seed(20)
p <- 10
mu <- rep(0, p)
X <- matrix(rnorm(p * (p-1)), p, p-1)
X <- cbind(X, X[, p-1])
S <- crossprod(X)
x <- rep(0, p)
err <- try(y <- dmvnorm(x, mu, S, log = TRUE))
cat(err)

options(op)
