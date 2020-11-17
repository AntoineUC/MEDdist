rmed=function(n,mu=rep(0,1),Sigma=diag(length(mu)),rho=0,nu=c(1,rep(0,length(mu)-1))){

  d=length(mu)
  
  normalizing_constant <- function(d = 2, # dimension, bivariate by default
                                   rho # denotes rho parameter: Gaussian distr with rho = 0
  ){
    nc <-
      if(d==1){
        (sqrt(1+rho)+sqrt(1-rho))/(2*sqrt(1-rho^2))
      }
    else if(d==2){
      1/(sqrt(1-rho^2))
    }
    else if(d==3){
      (sqrt(1+rho)-sqrt(1-rho))/(rho*sqrt(1-rho^2))
    }
    else if(d==4){
      2*(1-sqrt(1-rho^2))/(rho^2*sqrt(1-rho^2))
    }
    else
      stop("d is too large for my capacities...")
    nc*(2*pi)^(d/2)
    # nc <- log(beta(d/2-1/2,1/2))-d/2*log(2*pi)-log(simpadpt(
    #   function(psi_){sin(psi_)^(d-2)/(1+rho*cos(psi_))^(d/2)},0,pi))
    # nc <- 1/exp(nc)
  }
  
  upper_bound_constant <- function(d = 2, # dimension, bivariate by default
                                   rho # denotes rho parameter: Gaussian distr with rho = 0
  ){
    (2*pi/(1-rho))^(d/2)/normalizing_constant(d, rho)
  }
  
  dd <- 1:4
  rr <- (0:9)/10
  
  # pdf(file = "figures/acc-proba.pdf", width = 5, height = 5)
  # plot(NULL, xlim = c(1, max(dd)), ylim = c(0,1),
  #      xlab = expression(italic(d)), ylab = "Acceptance probability")
  # for(rho in rr){
  #   acc <- 1/upper_bound_constant(dd)
  #   lines(dd, acc, col =  rgb(rho,0,1-rho,0.8), lwd = .5)
  #   points(dd, acc, col =  rgb(rho,0,1-rho), pch = 19, lwd = .5)
  # }
  # dev.off()
  
  # density function for Sigma = I_d, mu = 0 and nu = e_1
  dexpectile <- function(x,
                         rho = 0 # denotes rho parameter: Gaussian distr with rho = 0
  ){
    d <- length(x)
    nc <- normalizing_constant(d, rho)
    norm_x_square <- sum(x^2)
    if(norm_x_square==0) 1/nc
    else{
      dot_prod <- x[1]/sqrt(norm_x_square)
      1/nc*exp(-norm_x_square/2*(1+rho*dot_prod))
    }
  }
  
  n_samp <- floor(n*upper_bound_constant(d, rho)*1.4)
  X <- rnorm(n_samp*d, 0, sd = 1/sqrt(1-rho))
  X <- matrix(X, ncol = d)
  U <- runif(n_samp)
  f_X <- apply(X, 1, function(x){dexpectile(x, rho)})
  g_X <- apply(X, 1, function(x){prod(dnorm(x, 0, sd = 1/sqrt(1-rho)))})
  accept <- (U*(2*pi/(1-rho))^(d/2)*(1/normalizing_constant(d,rho))*g_X < f_X) 
  if (d==1) {
    X <- X[which(accept==TRUE)]
    X <- X[1:(min(n,dim(X)[1]))]
    X
  } else{
    X <- X[which(accept==TRUE),]
    X <- X[1:(min(n,dim(X)[1])),]
    return(X)
  }

  
}