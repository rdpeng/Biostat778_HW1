dmvnorm <-
function (x, mu, S, log = TRUE) 
{
    l = tryCatch(chol(S), error = function(e) {
        message("S is not positive definite")
    })
    logdS = 2 * sum(log(diag(l)))
    n = nrow(x)
    Mx = x - rep(1, n) %*% t(mu)
    Y = forwardsolve(t(l), t(Mx))
    Mx3 = diag(crossprod(Y))
    k = length(mu)
    fx = -(k/2) * log(2 * pi) - (1/2) * logdS - (1/2) * Mx3
    if (log == FALSE) 
        fx = exp(fx)
    else fx = fx
    return(fx)
}
