
fastlm <- function(X, y, na.rm = FALSE) {
if (na.rm==TRUE){
Xrm<-which(is.na(X),arr.ind=TRUE)
yrm<-which(is.na(y))
getrid<-sort(unique(c(Xrm[,1],yrm)))
X<-X[-getrid,]
y<-y[-getrid]}
first<-as.matrix(solve(t(X)%*%X))
coefs<-as.matrix((first)%*%t(X)%*%y)
resids<-(diag(nrow(X))-(X%*%first%*%t(X)))%*%y
varyest<-(t(resids)%*%resids)/(nrow(X)-ncol(X))
mylist<-list("coefficients"=coefs,"vcov"=first*as.numeric(varyest))
return(mylist)
}