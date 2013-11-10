#define function
fastlm<-function(X, y, na.rm = FALSE) {
        n<-length(y)
        p<-ncol(X)
    #check argument na.rm
        if (na.rm==TRUE){
                Z=cbind(X,y)
                X=X[complete.cases(Z),]
                y=as.matrix(y[complete.cases(Z)])
        }
    
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
        cov_beta<-chol2inv(Q)*as.numeric(crossprod(y-X%*%betahat,y))/(n-p)
    
        return(list(coeffients=betahat,vcov=cov_beta))
}