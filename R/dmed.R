dmed=function(x,mu=rep(0,length(rbind(x)[,1])),Sigma=diag(length(rbind(x)[,1])),rho=0,nu=rep(1,length(rbind(x)[,1]))/sqrt(length(rbind(x)[,1]))){
  
  
  if(rho>=1 || rho<0){
    stop("rho must be between 0 and 1.")
  }
  
  if(det(Sigma)==0){
    stop("Sigma is singular.")
  }
  
  if(length(mu)!=length(rbind(x)[,1])){
    stop("x and mu must have the same length.")
  }
  
  if(length(nu)!=length(rbind(x)[,1])){
    stop("x and nu must have the same length.")
  }
  
  if(length(mu)!=length(rbind(x)[,1])){
    stop("x and mu must have the same length.")
  }
  
  #if(is.matrix(Sigma)!=F){
  #  stop("Sigma must be a matrix.")
  #}
  
  if(sqrt(length(Sigma))!=length(rbind(x)[,1])){
    stop("x and Sigma must have the same length.")
  }
  
  return(det(Sigma)^(-1/2)*(2*pi)^(-length(rbind(x)[,1])/2)*2^(1-length(rbind(x)[,1]))*(sqrt(1-rho)+sqrt(1+rho))^(length(rbind(x)[,1]))/(1+1/sqrt(1-rho^2))*exp(-0.5*diag(t(x-matrix(rep(mu,length(rbind(x)[1,])),length(rbind(x)[,1]),length(rbind(x)[1,])))%*%solve(Sigma)%*%(x-matrix(rep(mu,length(rbind(x)[1,])),length(rbind(x)[,1]),length(rbind(x)[1,]))))-0.5*rho*diag(t(x-matrix(rep(mu,length(rbind(x)[1,])),length(rbind(x)[,1]),length(rbind(x)[1,])))%*%solve(Sigma)%*%matrix(rep(nu,length(rbind(x)[1,])),length(rbind(x)[,1]),length(rbind(x)[1,])))*sqrt(diag(t(x-matrix(rep(mu,length(rbind(x)[1,])),length(rbind(x)[,1]),length(rbind(x)[1,])))%*%solve(Sigma)%*%(x-matrix(rep(mu,length(rbind(x)[1,])),length(rbind(x)[,1]),length(rbind(x)[1,])))))))
}
