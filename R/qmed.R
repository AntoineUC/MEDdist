qmed=function(probs,mu=0,Sigma=1,rho=0,nu=1){
  
  if(max(probs)>=1 || min(probs)<=0){
    stop("probs must be between 0 and 1.")
  }
  
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
  
  if(length(Sigma)>1){
    stop("Sigma must be a real.")
  }
  
  find_root=function(tau){
    
    fpar=function(x){
      return(pmed(x,mu=mu,Sigma=Sigma,rho=rho,nu=nu)-tau)
    }
    
    return(uniroot(fpar,c(-40,40))$root)
    
  }
  
  return(mapply(find_root,probs))
  
}
