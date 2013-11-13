#define function
fastlm<-function(X, y, na.rm=FALSE) {        
    #check argument na.rm
        if (na.rm!=FALSE) {
                Z=cbind(X,y)
                X=X[complete.cases(Z),]
                y=as.matrix(y)[complete.cases(Z)]
        }
        n<-length(y)
        p<-ncol(X)
    #calculating transpose(X)%*%X
        A<-crossprod(X)
    #calculating transpose(X)%*%y    
        C<-crossprod(X,y)
    
    #cholesky decomposition
        Q<-chol(A)
   
    #solve betahat
        temp1<-forwardsolve(t(Q),C) 
        betahat<-backsolve(Q,temp1) 
    
    #calculate covirance of beta
    #note that t(e)%*%e=t(e)%*%y=t(y)%*%y-t(y)%*%X%*%betahat
    #the second and the third expression is almost the same in my computer
        cov_beta<-chol2inv(Q)*as.numeric(crossprod(y)-crossprod(X%*%betahat))/(n-p)
    
        return(list(coefficients=betahat,vcov=cov_beta))
}