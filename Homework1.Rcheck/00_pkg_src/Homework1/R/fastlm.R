fastlm<-function(X, y, na.rm=FALSE)
{
      ## calculate estimated beta
      b<-crossprod(X,y)
      XX<-crossprod(X)
      cx<-chol(XX)
      temp<-forwardsolve(t(cx),b)
      beta<-backsolve(cx,temp)
      
      ## calculate covariance matrix of estimated beta
      n<-length(y)
      p<-ncol(X)
      e<-(y-tcrossprod(X,t(beta)))
      e_var<-sum(e^2)/(n-p)
      vcov<-e_var*chol2inv(cx)
      
      ## return result list
      result<-list(coefficients=beta,vcov=vcov)
}

