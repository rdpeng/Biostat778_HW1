fastlm <-
function(X, y, na.rm = FALSE) {
      if (na.rm == TRUE)
      {
        cb = cbind(x,y)
        X = X[complete.cases(cb),]
        y = as.matrix(y[complete.cases(cb)])
      }
     
        n=nrow(X)
        p=ncol(X)
        l=chol(crossprod(X,X))
        cp=crossprod(X,y)
        cof=forwardsolve(t(l),forwardsolve(t(l),cp),transp=TRUE)
        vcov<-as.numeric((crossprod(y,y)-crossprod(cof,cp)))*chol2inv(l)/(n-p)
        list(cof,vcov)## Your code goes here
}
