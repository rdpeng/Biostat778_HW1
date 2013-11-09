fastlm<-function(X, y, na.rm=FALSE)
{
      if(na.rm==T)
      {
            notna<-complete.cases(X,y)
            X<-X[notna,]
            y<-y[notna]
      }
      
      ## calculate estimated beta
      b<-crossprod(X,y)
      cx<-chol(crossprod(X))
      temp<-forwardsolve(t(cx),b)
      beta<-as.vector(backsolve(cx,temp))
      
      ## calculate covariance matrix of estimated beta
      e_var<-(sum(y^2)-sum(b*beta))/(length(y)-length(beta))
      vcov<-e_var*chol2inv(cx)
      
      ## return result list
      list(coefficients=beta,vcov=vcov)
}

