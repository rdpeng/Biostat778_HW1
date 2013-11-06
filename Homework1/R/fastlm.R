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
    y <- as.vector(y)
    if (na.rm)
    {
        to_keep <- complete.cases(X) & complete.cases(y)
        X <- X[to_keep,]
        y <- y[to_keep]
    }
    
    xx <- crossprod(X)
    xy <- crossprod(X,y)
    root <- chol(xx)
    rootbeta <- forwardsolve(t(root),xy)
    beta_hat <- backsolve(root, rootbeta)
    
    k <- ncol(root)
    rooti <- backsolve(root, diag(k))
    xx_inv <- tcrossprod(rooti)
    beta_hat_vcov <- xx_inv * var(y)

    list(coefficients=beta_hat, vcov=beta_hat_vcov)
}