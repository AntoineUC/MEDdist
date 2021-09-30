baymode=function(X,niter=100){

  n=length(X[,1])

  d=length(X[1,])

  med_dat <- list(n = n,
                  d = d,
                  Y = X
  )

  fit <- sampling(object = stanmodels$mode, data = med_dat,
              chains = 2, iter = niter)

  return(list(mu=extract(fit)$mu,muhat=apply(extract(fit)$mu,2,mean)))

}



