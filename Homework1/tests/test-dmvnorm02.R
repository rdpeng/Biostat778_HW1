## Check `log = FALSE'

library(Homework1)

mu <- rep(0, 10)
S <- diag(2, 10)
x <- matrix(rbind(rep(0, 10), rep(1, 10)), 2, 10)
y <- dmvnorm(x, mu, S, log = FALSE)
print(y, digits = 15)
y <- dmvnorm(x, mu, S, log = TRUE)
print(y, digits = 15)
