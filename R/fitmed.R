fitmed=function(x,mu0=colMeans(rbind(x)),Sigma0=var(x),rho0=0.1,nu0=sample(c(1,rep(0,length(rbind(x)[1,])-1))),maxit=100){
  
  mu_ <- mu0
  sigma_ <- Sigma0
  rho_ <- rho0
  nu_ <- nu0
  d_=length(rbind(x)[1,])
  n_=length(rbind(x)[,1])
  x_=x
  
  log_like <- function(mu_,sigma_,rho_,nu_) {
    sigma_inv_ <- ginv(sigma_)
    sigma_sqrt_inv_ <- sqrtm(sigma_)$Binv
    ll_ <- 0
    ll_ <- ll_ + n_*log(beta(d_/2-1/2,1/2))
    ll_ <- ll_ - n_*d_/2*log(2*pi) - n_/2*log(det(sigma_))
    ll_ <- ll_ - n_*log(simpadpt(
      function(psi_){sin(psi_)^(d_-2)/(1+rho_*cos(psi_))^(d_/2)},0,pi))
    # ll_ <- ll_ + n_*log(quadl(
    #   function(tt) {
    #     (1-tt^2)^((d_-3)/2)*(1+rho_*tt)^(-d_/2)
    #   },
    #   -1,1
    # ))
    for (ii in 1:n_) {
      maha_ <- sqrt(t(x_[ii,]-mu_)%*%sigma_inv_%*%(x_[ii,]-mu_))
      ll_ <- ll_ - 1/2*maha_^2 - 
        1/2*rho_*maha_*t(x_[ii,]-mu_)%*%sigma_sqrt_inv_%*%nu_
    }
    return(ll_)
  }
  
  mu_obj <- function(para_) {
    mu_sub_ <- c(para_) 
    sigma_inv_ <- ginv(sigma_)
    sigma_sqrt_inv_ <- sqrtm(sigma_)$Binv
    obj_ <- 0
    for (ii in 1:n_) {
      maha_ <- sqrt(t(x_[ii,]-mu_sub_)%*%sigma_inv_%*%(x_[ii,]-mu_sub_))
      obj_ <- obj_ + maha_^2 +
        rho_*maha_*t(x_[ii,]-mu_sub_)%*%sigma_sqrt_inv_%*%nu_
    }
    return(obj_)
  }
  
  rho_obj <- function(para_) {
    rho_sub_ <- c(para_)
    sigma_inv_ <- ginv(sigma_)
    sigma_sqrt_inv_ <- sqrtm(sigma_)$Binv
    obj_ <- 0
    obj_ <- obj_ + 2*n_*log(simpadpt(
      function(psi_){sin(psi_)^(d_-2)/(1+rho_sub_*cos(psi_))^(d_/2)},0,pi))
    # obj_ <- obj_ + 2*n_*log(quadl(
    #   function(tt) {
    #     (1-tt^2)^((d_-3)/2)*(1+rho_sub_*tt)^(-d_/2)
    #   },
    #   -1,1
    # ))
    for (ii in 1:n_) {
      maha_ <- sqrt(t(x_[ii,]-mu_)%*%sigma_inv_%*%(x_[ii,]-mu_))
      obj_ <- obj_ + maha_^2 +
        rho_sub_*maha_*t(x_[ii,]-mu_)%*%sigma_sqrt_inv_%*%nu_
    }
    return(obj_)
  }
  
  nu_obj <- function(para_) {
    nu_sub_ <- c(para_) 
    sigma_inv_ <- ginv(sigma_)
    sigma_sqrt_inv_ <- sqrtm(sigma_)$Binv
    obj_ <- 0
    for (ii in 1:n_) {
      maha_ <- sqrt(t(x_[ii,]-mu_)%*%sigma_inv_%*%(x_[ii,]-mu_))
      obj_ <- obj_ + rho_*maha_*t(x_[ii,]-mu_)%*%sigma_sqrt_inv_%*%nu_sub_
    }
    return(obj_)
  }
  
  sigma_obj <- function(para_) {
    sigma_sub_ <- matrix(para_,d_,d_)
    sigma_inv_ <- ginv(sigma_sub_)
    sigma_sqrt_inv_ <- sqrtm(sigma_sub_)$Binv
    obj_ <- 0
    obj_ <- obj_ + n_*log(det(sigma_sub_))
    for (ii in 1:n_) {
      maha_ <- sqrt(t(x_[ii,]-mu_)%*%sigma_inv_%*%(x_[ii,]-mu_))
      obj_ <- obj_ + maha_^2 +
        rho_*maha_*t(x_[ii,]-mu_)%*%sigma_sqrt_inv_%*%nu_
    }
    return(obj_)
  }
  
  
  mani_def_nu_ <- get.sphere.defn(d_)
  mani_params_nu_ <- get.manifold.params(IsCheckParams = FALSE)
  solver_params <- get.solver.params(Tolerance = 1e-4, Max_Iteration = 10)
  
  deriv_params_nu_ <- get.deriv.params()
  
  
  mod <- Module("ManifoldOptim_module", PACKAGE = "ManifoldOptim")
  prob_nu_ <- new(mod$RProblem, nu_obj)
  
  mani_def_sigma_ <- get.spd.defn(d_)
  mani_params_sigma_ <- get.manifold.params(IsCheckParams = FALSE)
  
  deriv_params_sigma_ <- get.deriv.params()
  
  prob_sigma_ <- new(mod$RProblem, sigma_obj)
  
  for (rr in 1:maxit) {
    
    opt_mu_ <- optim(mu_,mu_obj)
    mu_ <- opt_mu_$par
    # ll_ <- log_like(mu_,sigma_,rho_,nu_)
    # print(ll_)
    
    opt_sigma_ <- manifold.optim(prob_sigma_, mani_def_sigma_, c(sigma_), method = "LRTRSR1",
                                 mani.params = mani_params_sigma_, solver.params = solver_params, deriv.params = deriv_params_sigma_)
    sigma_ <- matrix(opt_sigma_$xopt,d_,d_)
    # ll_ <- log_like(mu_,sigma_,rho_,nu_)
    # print(ll_)
    
    opt_rho_ <- optim(rho_,rho_obj,method = 'Brent',lower = 0, upper = 1-1e-6)
    rho_ <- opt_rho_$par
    # ll_ <- log_like(mu_,sigma_,rho_,nu_)
    # print(ll_)
    
    opt_nu_ <- manifold.optim(prob_nu_, mani_def_nu_, nu_, method = "LRTRSR1",
                              mani.params = mani_params_nu_, solver.params = solver_params, deriv.params = deriv_params_nu_)
    nu_ <- c(opt_nu_$xopt)
    
    ll_ <- log_like(mu_,sigma_,rho_,nu_)
    
  } 
  paramlist=list("mu"=mu_,"Sigma"=sigma_,"rho"=rho_,"nu"=nu_)
  
  return(paramlist)
  
}