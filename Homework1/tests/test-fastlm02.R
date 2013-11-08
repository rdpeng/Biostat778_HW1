## Remove missing values

library(Homework1)
library(datasets)
data(airquality)

X <- data.matrix(airquality[, -1])
y <- airquality$Ozone

fit <- fastlm(X, y, na.rm = TRUE)
print(drop(fit$coefficients))
print(fit$vcov)
