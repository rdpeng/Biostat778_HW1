#define function
fastlm<-function(X, y, na.rm = FALSE) {
        n<-length(y)
        p<-ncol(X)
    #check na.rm
        if (na.rm==TRUE){
                Z=cbind(X,y)
                X=X[complete.cases(Z),]
                y=as.matrix(y[complete.cases(Z)])
        }
    
    #calculating transpose(X)%*%X
        A<-crossprod(X)
        C<-crossprod(X,y)
    
    #cholesky decomposition
        Q<-chol(A)
   
    #solve betahat
        temp1<-forwardsolve(t(Q),C) 
        betahat<-backsolve(Q,temp1) 
    
    #calculate covirance of beta
        cov_beta<-chol2inv(Q)*as.numeric(crossprod(y-X%*%betahat)/(n-p))
    
        return(list(coeffients=betahat,vcov=cov_beta))
}