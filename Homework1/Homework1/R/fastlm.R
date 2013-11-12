#' Fast linear regression
#'
#' @param X an \code{n x p} design matrix
#' @param y a numeric vector of length \code{n}
#' @param na.rm a boolean to indicate if \code{NA} values should be removed
#' @return \code{fastlm()} returns a list with the following components: \item{coefficients}{A numeric vector of length \code{p} containing the regression coefficients estimated by maximum likelihood.} \item{vcov}{the \code{p x p} covariance matrix of the estimated regression coefficients.}
#' @export
#' @examples
#' X <- matrix(rnorm(10),ncol=2)
#' y <- rnorm(5)
#' fastlm(X, y)
fastlm <- function(X, y, na.rm=FALSE) 
{   
    if (na.rm)
    {
        to_keep <- complete.cases(X,y)
        X <- X[to_keep,]
        y <- y[to_keep]
    }
    
    n <- nrow(X)
    p <- ncol(X)
    
    xx <- crossprod(X)
    xy <- crossprod(X,y)
    root <- chol(xx)
    rootbeta <- forwardsolve(t(root),xy)
    beta_hat <- backsolve(root, rootbeta)
    
    root_inv <- backsolve(root, diag(p))
    xx_inv <- tcrossprod(root_inv)
    s2 <- (sum(y*y) - sum(xy*beta_hat)) / (n - p)
    beta_hat_vcov <- xx_inv * s2
    
    list(coefficients=beta_hat, vcov=beta_hat_vcov)
}