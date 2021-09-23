pmed=function(x,mu=0,Sigma=1,rho=0,nu=1){
  
  
  if(rho>=1 || rho<0){
    stop("rho must be between 0 and 1.")
  }
  
  if(Sigma<=0){
    stop("Sigma must be positive.")
  }
  
  if(length(mu)>1){
    stop("mu must be a real.")
  }
  
  if(length(nu)>1){
    stop("nu must be a real.")
  }
  
  if(length(x)>1){
    stop("x must be a real.")
  }
  
  #if(is.matrix(Sigma)!=F){
  #  stop("Sigma must be a matrix.")
  #}
  
  if(length(Sigma)>1){
    stop("Sigma must be a real.")
  }
  
  dmed1=function(x,mu=0,Sigma=1,rho=0,nu=1){
    cste=sqrt(1-rho^2)*((1+sqrt(1-rho^2))/2)^(-1/2)
    return(cste/(sqrt(2*pi)*sqrt(Sigma))*exp(-0.5*abs(x)/Sigma*(abs(x)+rho*x*nu/sqrt(Sigma))))
  }
  
  fpar=function(x){
    return(dmed1(x,mu=mu,Sigma=Sigma,rho=rho,nu=nu))
  }
  
  return(integrate(f=fpar,lower=-Inf,upper=x)$value)
  
}
