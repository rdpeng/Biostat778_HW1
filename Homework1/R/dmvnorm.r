

dmvnorm <- function(x, mu, S, log = TRUE) {
checkem<-eigen(S)$values
if(sum(checkem>0)!=length(checkem)){
stop("S is not positive definite")}
k<-length(x)
expterm<-((-1/2)*sqrt(mahalanobis(x,mu,S)))
dens<-exp(((-k/2)*log(2*pi))-((1/2)*log(det(S)))-expterm)
if(log==TRUE){
return(log(dens))}
return(dens)
}