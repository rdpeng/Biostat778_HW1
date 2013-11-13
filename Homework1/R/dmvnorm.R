#' Evaluate the density of a multivariate normal distribution
#'
#' @param x an \code{n x k} matrix to represent \code{n} vectors of length \code{k}
#' @param mu the mean of the distribution: a numeric vector of length \code{k}
#' @param S the variance-covariance matrix with dimensions \code{k x k}
#' @param log a boolean that indicates if the log of the density should be returned
#' @return numeric vector with length \code{n}, corresponding to the density of the \code{n} vectors in \code{x} evaluated under a multivariate normal distribution defined by mean \code{mu} and variance-covariance \code{S}
#' @export
#' @examples
#' x <- matrix(rnorm(10*9), ncol=9)
#' mu <- rep(0,9)
#' xg <- seq(0, 1, length = 3)
#' yg <- xg
#' g <- data.matrix(expand.grid(xg, yg))
#' D <- as.matrix(dist(g))
#' S <- exp(D * -1)
#' dmvnorm(x, mu, S)
dmvnorm <- function(x, mu, S, log=TRUE)
{
    k <- ncol(S)
    # find inverse of chol(Sigma)
    # test for Positive definiteness
    root <- tryCatch(
    {
        chol(S)
    },
    error=function(cond)
    {
        stop("S is not positive definite")
    })
    
    root_inv <- backsolve(root, diag(k))
    
    n_row <- nrow(x)
    z_vec <- rep(0, n_row)
    for (j in 1:n_row)
    {
        z <- crossprod(root_inv, (x[j,]-mu))
        z_vec[j] <- -0.5*as.numeric(crossprod(z))
    }
    
    log_val <- z_vec - k/2*base::log(2*pi) + sum(base::log(diag(root_inv)))
    if(log)
    {
        return(log_val)
    }
    else
    {
        return(exp(log_val))
    }
}