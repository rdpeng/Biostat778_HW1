## Check `log = FALSE'

library(Homework1)

mu <- rep(0, 100)
S <- diag(1, 100)
x <- rep(0, 100)
y <- dmvnorm(x, mu, S, log = TRUE)
print(y, digits = 15)
