bayexpectile=function(X,dir=c(1,rep(0,length(X[1,])-1)),rho=0.5,niter=100){
  
  setwd(.libPaths())
  
  setwd("MEDdist")

  if(rho>=1 || rho<0){
    stop("rho must be between 0 and 1.")
  }

  if(sum(dir^2)!=1){
    stop("dir must be on the unit circle.")
  }

  nu_tilde=dir

  n=length(X[,1])

  d=length(X[1,])

  med_dat <- list(n = n,
                  d = d,
                  rho=rho,
                  nu_tilde=nu_tilde,
                  Y = X
  )

  fit <- stan(file = 'expectile.stan', data = med_dat,
              chains = 2, iter = niter)

  return(list(mu=extract(fit)$mu,muhat=apply(extract(fit)$mu,2,mean)))

}



