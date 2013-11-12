fastlm <-
function (x, y, na.rm = FALSE) 
{
    if (na.rm == TRUE) {
        r = cbind(x, y)
        x = x[complete.cases(r), ]
        y = as.matrix(y[complete.cases(r)])
    }
    A = crossprod(x)
    B = crossprod(x, y)
    l = chol(A)
    beta = backsolve(l, forwardsolve(t(l), B))
    k = ncol(x)
    n = nrow(y)
    sigmaS = (crossprod(y, y) - crossprod(beta, B))/(n - k)
    covBeta = as.numeric(sigmaS) * chol2inv(l)
    list(beta, covBeta)
}
