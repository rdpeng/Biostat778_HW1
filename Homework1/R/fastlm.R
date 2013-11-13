fastlm<- function(X, y, na.rm = FALSE)
{
     
      if(is.matrix(X)==F) X<-as.matrix(X)
      #if(is.vector(y)==F) y<-as.vector(y)
      
      if(na.rm==T){
            ok<-complete.cases(X,y)
            X<-X[ok,]
            y<-y[ok]
      }
      
## compute QR-decomposition of X
qx <- qr(X)
## compute (X'X)^(-1) X'y
coef <- solve.qr(qx, y)
## degrees of freedom and standard deviation of residuals
df <- nrow(X)-ncol(X)
sigma2 <- sum((y - X%*%coef)^2)/df
## compute sigma^2 * (X'X)^-1
vcov <- sigma2 * chol2inv(qx$qr)
list(coefficients = coef,
vcov = vcov)
}
