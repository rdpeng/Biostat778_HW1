##HW1 Fast Linear Regression
fastlm=function(X,y,na.rm=FALSE){
        if(na.rm==TRUE){
                r=cbind(X,y)
                X=X[complete.cases(r),]
                y=as.matrix(y[complete.cases(r)])
        }
        
        ##Cholesky factorization for coefficients
        A=crossprod(X)
        B=crossprod(X,y)
        R=chol(A)
        Rbeta=forwardsolve(t(R),B)
        coefficients=backsolve(R,Rbeta)
        
        ##Calculate VCOV
        n=length(y)
        p=ncol(X)
        sigmahat2=crossprod(y-X%*%coefficients)/(n-p)
        Ainv=chol2inv(R)
        vcov=as.numeric(sigmahat2)*Ainv
        
        return(list(coefficients=coefficients,vcov=vcov))
}



