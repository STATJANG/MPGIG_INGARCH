#Necessary libraries.
library(GeneralizedHyperbolic)
library(gamlss)
library(tscount)
library(GIGrvg)
library(DistributionUtils)  







##Data generation with the parameterization mu=1.
old_mpgig_ingarch_sim = function(n=300,model,phi,alpha,d,B,A=NULL,burns=500){
  
  #set.seed(seedv)
  p = dim(B)[1]
  max.ab = max(model[["past_obs"]])
  nu = matrix(0,ncol = p,nrow = n+burns)
  dat = matrix(0,ncol = p,nrow = n+burns)
  gig_z=0
  b.ratio=besselRatio(phi,alpha,1)
  
  if(is.null(A)==F){
    for( i in 1:(n+burns) ){
      
      #gig_z[i] = gamlss.dist::rGIG(1,mu=1,sigma=1/sqrt(phi),nu=alpha)
      #gig_z[i] = GIGrvg::rgig(1,chi=phi/b.ratio,psi=phi*b.ratio,lambda=alpha) #mu=1
      gig_z[i] = GIGrvg::rgig(1,chi=phi,psi=phi,lambda=alpha) #mu=R(phi,alpha)
      if(i<=max.ab) {
        dat[i,] = rpois(p,exp(nu[i,])*gig_z[i])
      }else{
        nu[i,] = d + A%*%c(nu[i-model[["past_mean"]],]) + B%*%log(c(dat[i-model[["past_obs"]],])+1) 
        dat[i,] = rpois(p,exp(nu[i,])*gig_z[i])
      }
    } # end loop over i
  }else{
    
    for( i in 1:(n+burns) ){
      #gig_z[i] = gamlss.dist::rGIG(1,mu=1,sigma=1/sqrt(phi),nu=alpha)
      #gig_z[i] = GIGrvg::rgig(1,chi=phi/b.ratio,psi=phi*b.ratio,lambda=alpha) #mu=1
      gig_z[i] = GIGrvg::rgig(1,chi=phi,psi=phi,lambda=alpha) #mu=R(phi,alpha)     
      if(i<=max.ab) {
        dat[i,] = rpois(p,exp(nu[i,])*gig_z[i])
      }else{
        nu[i,] = d +B%*%log(c(dat[i-model[["past_obs"]],])+1) 
        dat[i,] = rpois(p,exp(nu[i,])*gig_z[i])
      }
    }
    
  }
  return( dat[-(1:burns),] )
}


# Approximation of Bessel in case it diverges due to a large argument(phi).
approx_logK=function(x,nu){
  log(pi/2/nu)/2-nu*log(exp(1)*x/2/nu)
}




#parameterizatino mu=R(phi,alpha)
old_INGARCH_init = function(dat,model=list(past_mean=1,past_obs=1)){
  
  dat_n = dim(dat)[1] #the number of observation
  dat_p = dim(dat)[2] #the dimension of multivariate time series; 
  init = new_INGARCH_init(dat,model)
  init_phi = init_alpha = rep(NA,dat_p) 
  
  if(is.null(model[["past_mean"]])==T){ #in'A'rch
    
    for(p in 1:dat_p){
      tryCatch(
        expr = {
          init_GIG = gamlss(dat[(3:dat_n),p]~1+log(1+dat[2:(dat_n-1),p]),family="SICHEL",control=gamlss.control(trace=FALSE))
          init_phi[p] = max(exp(-init_GIG$sigma.coefficients),0.1)
          init_alpha[p] = init_GIG$nu.coefficients
          
        },
        error = function(e){ 
        }
      )
    }
    tmp_est = c( median(init_phi,na.rm = T),median(init_alpha,na.rm = T) )
    tmp_d = init$estimate[3:(3+dat_p)]/besselRatio(tmp_est[1],tmp_est[2],1)
    #tmp_d = init$estimate[3:(3+dat_p)]/mean(GIGrvg::rgig(10000,tmp_est[1],tmp_est[1],lambda=tmp_est[2]))
    tmp_est = c(tmp_est,tmp_d,init$estimate[-(1:(3+dat_p))])
    z=list(estimate=tmp_est,fitted = init$fitted, residuals = init$resid)
    class(z) = "ingarch_init"
    return(z)
  }else{
    
    for(p in 1:dat_p){
      tryCatch(
        expr = {
          init_GIG = gamlss(dat[(3:dat_n),p]~1+log(1+dat[2:(dat_n-1),p]),family="SICHEL",control=gamlss.control(trace=FALSE))
          init_phi[p] = max(exp(-init_GIG$sigma.coefficients),0.1)
          init_alpha[p] = init_GIG$nu.coefficients
        },
        error = function(e){ 
        }
      )
    }
    tmp_est = c( median(init_phi,na.rm = T),median(init_alpha,na.rm = T) )
    tmp_d = init$estimate[3:(3+dat_p)]/besselRatio(tmp_est[1],tmp_est[2],1)
    tmp_est = c(tmp_est,tmp_d,init$estimate[-(1:(3+dat_p))])
    z=list(estimate=tmp_est,fitted = init$fitted, residuals = init$resid)
    class(z) = "ingarch_init"
    return(z)
  }
}

#Obtain initial values (mu=1 parameterization)
new_INGARCH_init = function(dat,model=list(past_mean=1,past_obs=1),distr = "MSICHEL"){#'dat' here is the full data.
  
  
  dat_n = dim(dat)[1] #the number of observation
  dat_p = dim(dat)[2] #the dimension of multivariate time series; 
  bb = length(model[["past_obs"]])
  max.ab = max(model[["past_obs"]])
  
  ##Initial values for d, A, B.
  init_d = rep(0,dat_p) 
  init_B = matrix(0,ncol=dat_p*bb ,nrow=dat_p)
  init_lam = init_resid=matrix(0,ncol=dat_p,nrow=dat_n)
  disp = 0
  
  if(is.null(model[["past_mean"]])==T){ #in'A'rch
    #num_par = 2+dat_p+ (aa+bb)*dat_p^2 # the number of parameters.
    
    for(p in 1:dat_p){
      ##Initial value for ARCH coefficient, 
      init_val = tsglm(ts = dat[,p], model = model,link = "log",distr = "nbinom")
      init_d[p] = unname(init_val$coefficients[1])
      init_B[p,(1:bb)+(p-1)*bb] = unname(init_val$coefficients[(1:bb)+1])
      disp[p] = unname(init_val$sigmasq)
      init_lam[,p] = init_val$fitted.values
      init_resid[,p] = init_val$residuals
    }
    z=list(estimate=c(mean(disp),-0.5,init_d,init_B),
           fitted = init_lam[-c(1:max.ab),], residuals = init_resid[-c(1:max.ab),])
    if(distr=="MNB"){z=list(estimate=c(1/mean(disp),init_d,init_B),
                            fitted = init_lam[-c(1:max.ab),], residuals = init_resid[-c(1:max.ab),])}
    class(z) = "ingarch_init"
    return(z)
  }else{
    
    aa = length(model[["past_mean"]])
    init_A = matrix(0,ncol=dat_p*aa,nrow=dat_p)
    for(p in 1:dat_p){
      ##Initial value for GARCH coefficient, 
      init_val = tsglm(ts = dat[,p], model = model,link = "log",distr = "nbinom")
      init_d[p] = unname(init_val$coefficients[1])
      init_A[p,(1:aa)+(p-1)*aa] = unname(init_val$coefficients[(1:aa)+1+bb])
      init_B[p,(1:bb)+(p-1)*bb] = unname(init_val$coefficients[(1:bb)+1])
      disp[p] = unname(init_val$sigmasq)
      init_lam[,p] = init_val$fitted.values
      init_resid[,p] = init_val$residuals
    }
    z=list(estimate=c(mean(disp),-0.5,init_d,init_B,init_A),
           fitted = init_lam[-c(1:max.ab),], residuals = init_resid[-c(1:max.ab),])
    if(distr=="MNB"){z=list(estimate=c(1/mean(disp),init_d,init_B,init_A),
                            fitted = init_lam[-c(1:max.ab),], residuals = init_resid[-c(1:max.ab),])}
    if(distr=="Poisson"){z=list(estimate=c(init_d,init_B,init_A),
                                fitted = init_lam[-c(1:max.ab),], residuals = init_resid[-c(1:max.ab),])}
    class(z) = "ingarch_init"
    return(z)
  }
}
print.ingarch_init = function(x){print(x["estimate"])}

    






#Compute the mean process of log-linear model.
#(regardless of mu parameterization.)
mean_process = function(dat,pars,model=list(past_mean=1,past_obs=1),Xcov=NULL){
  #'dat' here is the full data.
  #dat should be at least two-dim matrix.
  #'model' should be a format like 'model = list(past_mean=c(1,12),past_obs=c(1,12))'
  #Xcov should be a matrix even if one-dim.
  #pars = c(d,B,A) or c(d,B,A,C)
  #All the matrix is column-wise when converted into a vector.
  
  #Errors
  if(is.null(model[["past_obs"]])==T){stop("At least one observation term must be included.")}
  
  n=dim(dat)[1]
  p=dim(dat)[2]
  bb  = length(model[["past_obs"]])
  max.ab = max(model[["past_obs"]])
  nu = matrix(0,ncol=p,nrow=n) #container
  d = c(pars[1:p])
  B = matrix(pars[(1:(bb*p^2))+p],nrow=p,ncol=p*bb)
  
  if(is.null(model[["past_mean"]])==F){   #########in'G'arch model
    if(all(model[["past_mean"]]%in%model[["past_obs"]])==F){stop("Any past mean term cannot be used independently.")}
    
    aa = length(model[["past_mean"]])  
    A = matrix(pars[(1:(aa*p^2))+p+bb*p^2],nrow=p,ncol=p*aa)
    
    if(is.null(Xcov)==F){ #With covariate.
      px = dim(Xreg)
      C = matrix(pars[-(1:((aa+bb)*p^2+p) )],ncol=px,nrow=p)  
      for(t in (max.ab+1):n){
        nu[t,] = d + A%*%c(nu[t-model[["past_mean"]],1:p]) + B%*%log(as.vector(dat[t-model[["past_obs"]],1:p])+1) + C%*%Xcov[t,]
      }
    }else{ #Without Covariate
      for(t in (max.ab+1):n){
        nu[t,] = d + A%*%c(nu[t-model[["past_mean"]],1:p]) + B%*%log(as.vector(dat[t-model[["past_obs"]],1:p])+1)
      }
    }
    
    
  }else{  ############in'A'rch model.
    
    
    
    if(is.null(Xcov)==F){ #With covariate.
      px = dim(Xreg)
      C = matrix(pars[-(1:((aa+bb)*p^2+p) )],ncol=px,nrow=p)  
      for(t in (max.ab+1):n){
        nu[t,] = d + B%*%log(as.vector(dat[t-model[["past_obs"]],1:p])+1)+ C%*%Xcov[t,]
      }
    }else{ #Without Covariate
      for(t in (max.ab+1):n){
        nu[t,] = d + B%*%log(as.vector(dat[t-model[["past_obs"]],1:p])+1)
      }
    }
    
  }
  return(exp(nu[-(1:max.ab),]))
}



#Pmf evaluation of multivariate Sichel distr for a single observation
#(parameterization;  mu=R(phi,alpha))
old_mSichel_pairs = function(y2, phi2, alpha2, lambda2){
  
  term1 = log(besselK( x=sqrt(phi2*2*sum(lambda2)+phi2^2  ) , nu=alpha2 +sum(y2)))
  #if(term1==Inf){ term1 =  log(0.5) + lgamma(sum(y2)+alpha2+1) - 
  #  (sum(y2)+alpha2)*log(sqrt(phi2*2*sum(lambda2)+phi2^2)/2)}
  
  #this formula is accurate.
  if(is.infinite(term1)){ term1 = approx_logK(sqrt(phi2*2*sum(lambda2)+phi2^2  ), nu=alpha2 +sum(y2))}
  
  
  term2 = log(besselK( x=phi2 , nu=alpha2 ))
  term3 = sum( y2*log(lambda2)-lgamma(y2+1) )
  term4 = ( (sum(y2)+alpha2) )*( log(phi2)- 0.5*log( 2*sum(lambda2)*phi2 +  phi2^2)) 
  return(term1-term2+term3+term4)
} 


#likelihood of MPGIG_INGARCH model with contemporaneous multivariate distr = multivariate Sichel.
#(parameterization;  mu=R(phi,alpha))
old_mSichel_f = function(parms,dat,model){
  
  lambdam = mean_process(dat,parms[-(1:2)],model) #dim(lambdam)[2] = n-max.ab
  dat_p = dim(dat)[2]
  max.ab = max(model[["past_obs"]])
  phi = parms[1] # \phi
  alpha = parms[2] # \alpha
  p = dim(dat)[2] # dimension of MV series
  
  #outp = 0; for(i in 2:n) outp = outp + mSichel_pairs( dat[i,] , phi, alpha , lambdam[,i-1] )
  outp = sum(apply(cbind(dat[-(1:max.ab),],lambdam),1,
                   function(x){old_mSichel_pairs(y2=x[1:dat_p],phi2=phi,alpha2=alpha,lambda2=x[-(1:dat_p)])}))
  
  return(-outp) # negative log likelihood
} 






##########Func's for EM
##########Func's for EM
##########Func's for EM

#Q2 is common irrelevant to the choice of a mixing distribution. (Poisson-based)
common_Q2 = function(dat,cond_z,pars,model){   
  #cond_z = a vector of E(Z_t|y_t).
  
  max.ab = max(model[["past_obs"]])
  lambdam = mean_process(dat,pars,model)
  
  result = -sum(rowSums(lambdam)*cond_z) + sum(rowSums(log(lambdam)*dat[-(1:max.ab),,drop=F]))
  return(-result)
}






#####EM under older parameterization###
# Monte Carlo approx of ( E(g(Z)|y) )
#(parameterization;  mu=R(phi,alpha))
old_GIG_cond_mc = function(tmp_dat, phi, alpha,lambda){
  #Based on the new parametrization GIG(mu=1,phi,alpha)
  tmp1 = log(besselK(phi,alpha))
  tmp2 = log(besselK(phi,alpha+1))
  cond = rowSums(cbind(is.infinite(tmp1),is.infinite(tmp2)))>0
  if(sum(cond)>0){
    tmp1[cond] = approx_logK(phi[cond],alpha[cond])
    tmp2[cond] = approx_logK(phi[cond],alpha[cond]+1)
  }
  
  eta = exp(tmp1)/exp(tmp2) #(1/R1)
  
  phi2 = sqrt(phi^2+2*phi*rowSums(lambda))
  alpha2 = alpha + rowSums(tmp_dat)
  eta2 = phi/phi2
  tmp_par = cbind(eta2*phi2,phi2/eta2,alpha2)  #(chi=b, psi=a, lambda=alpha)
  
  result = matrix(0,nrow=dim(tmp_dat)[1],ncol=3)
  for(i in 1:dim(tmp_dat)[1]){
    tmp = GIGrvg::rgig(10000,chi=tmp_par[i,1],psi=tmp_par[i,2],lambda=tmp_par[i,3])
    result[i,] = c( mean(tmp), mean(1/tmp), mean(log(tmp)) )
  }
  return(list(mean_vector = colMeans(result), cond_z = result[,1]))
} 


#(parameterization;  mu=R(phi,alpha))
old_GIG_Q1 = function(phi,alpha,cond_exp){  #cond_exp = c(E(Z|y),E(Z^-1|y),E(log(Z)|y))
  
  tmp = log(besselK(phi,alpha))
  if(is.nan(tmp)){tmp=approx_logK(phi,alpha)}
  if(is.infinite(tmp)){tmp=approx_logK(phi,alpha)}
  
  result = 
    alpha*cond_exp[3]-tmp-0.5*phi*sum(cond_exp[1:2])
  return(-result)
}



#(parameterization;  mu=R(phi,alpha))
old_INGARCH_EM = function(dat,init=NULL,model,tol_em = 0.001,max_iter=50, trace=F,method="BFGS"){
  
  init_val = init
  if( is.null(init) == T ){ init_val =old_INGARCH_init(dat,model)$estimate}
  
  dat_n = dim(dat)[1] #the number of observation
  dat_p = dim(dat)[2] #the dimension of multivariate time series; 
  
  bb = length(model[["past_obs"]])
  aa = length(model[["past_mean"]])
  max.ab = max(model[["past_obs"]])
  tmp_dat = dat[-c(1:max.ab),]
  
  #if( is.null(model[["past_mean"]])==F ){
  #  aa= length(model[["past_mean"]])
  #  num_par = 2+dat_p+(bb*dat_p^2 )+(aa*dat_p^2 )
  #}
  num_par = 2+dat_p + (bb*dat_p^2 ) + (aa*dat_p^2 )
  
  old.est = new.est = init_val; 
  dist=0;stored_loglik = 0;dist_em=10;iter=1;
  while(dist_em>tol_em){
    iter=iter+1
    #E-step
    lambda = mean_process(dat,old.est[-(1:2)],model)
    conds = old_GIG_cond_mc(tmp_dat,old.est[1],old.est[2],lambda)
    cond_z = conds$cond_z
    cond_exp = conds$mean_vector
    #cond_z = GIG_cond_mc_Z1(tmp_dat,old.est[1],old.est[2],lambda) #a vector of expected_Z
    #cond_exp =  GIG_cond_mc(tmp_dat,old.est[1],old.est[2],lambda)  #vector of three conditional moments
    #####The end of E-step#####
    
    #####M-step#####
    new.est[1:2] =  #Based on the parameterization GIG(mu=1,phi,alpha)
      optim(c(1/old.est[1],old.est[2]),
            function(pars){old_GIG_Q1(cond_exp,phi=1/pars[1],alpha=pars[2])},method=method)$par
    new.est[1] = 1/new.est[1]
    
    new.est[-c(1:2)] = 
      optim(old.est[-c(1:2)], function(pars){
        common_Q2(dat,cond_z,pars,model)},method=method)$par
    
    #####The end of M-step#####
    
    dist_em = max(abs(old.est-new.est))
    stored_loglik[iter] = old_mSichel_f(new.est,dat,model)
    if(trace==TRUE){
      if(iter==2){  stored_loglik[iter-1] = old_mSichel_f(old.est,dat,model) 
      cat("Initial value with loglik=",stored_loglik[iter-1],"\n")
      }
      cat("\n","iter=",iter,"\n",round(new.est,4),"\n","THE DIST IS ",dist_em,"\n","THE LOGLIK IS ",stored_loglik[iter],"\n")
      cat("dist=",dist_em,"and loglik=",stored_loglik[iter],"\n")
    }
    if(iter>max_iter ){break()}
    old.est = new.est
  }
  
  
  if( is.null(model[["past_mean"]])==F ){ #GARCH
    results = list(estimate_vector=new.est, log_lik = -stored_loglik[iter], iter=iter, conv_criterion = dist,fitted=lambda, 
                   residuals = tmp_dat-lambda*mean(cond_z), expected_Z = cond_z, AIC = stored_loglik[iter]+2*num_par,
                   BIC = stored_loglik[iter]+num_par*log(dim(tmp_dat)[1]),
                   estimate = list(phi=new.est[1],alpha=new.est[2],d=new.est[(1:dat_p)+2],B=matrix(new.est[(1:(bb*dat_p^2))+2+dat_p],nrow=dat_p),
                                   A=matrix(new.est[(1:(aa*dat_p^2))+2+dat_p+bb*dat_p^2],nrow=dat_p)),sample_size=dat_n, num_par = num_par ,model=model)
  }else{ #ARCH
    results = list(estimate_vector=new.est, log_lik = -stored_loglik[iter], iter=iter, conv_criterion = dist,fitted=lambda, 
                   residuals = tmp_dat-lambda*mean(cond_z), expected_Z = cond_z, AIC = stored_loglik[iter]+2*num_par,
                   BIC = stored_loglik[iter]+num_par*log(dim(tmp_dat)[1]),
                   estimate = list(phi=new.est[1],alpha=new.est[2],d=new.est[(1:dat_p)+2],B=matrix(new.est[(1:(bb*dat_p^2))+2+dat_p],nrow=dat_p),A=NULL),
                   sample_size=dat_n, num_par = num_par, model=model)  
  }    
  class(results) = "ingarch_em"
  return(results)
}
print.ingarch_em = function(x){print(x[c("estimate","log_lik","AIC","BIC","iter","conv_criterion")])}



##Parametric bootsrapping. - mu=R(phi,alpha)
old_ingarch_pb = function(ingarch_em,pb_n=500,burnin=500){
  if(class(ingarch_em)!="ingarch_em") stop("The class of input must be 'ingarch_em'.") 
  x = ingarch_em
  pb.mat = matrix(0,ncol=x$num_par,nrow=pb_n)
  pb.cor = 0
  pb.mean =pb.var= pb.kurtosis=pb.skew=pb.skew.nb = matrix(ncol=nrow(x$estimate$B),nrow=pb_n)
  pb.quantile1 = matrix(ncol=6,nrow=pb_n)
  pb.quantile2 = matrix(ncol=6,nrow=pb_n)  
  
  pb=1
  while(pb<=pb_n){
    tryCatch(
      expr = {
        pb.dat = old_mpgig_ingarch_sim(n=x$sample_size,model = x$model,phi = x$estimate$phi,alpha = x$estimate$alpha,d = x$estimate$d,B=x$estimate$B,A=x$estimate$A,burns=burnin)
        pb.em = old_INGARCH_EM(pb.dat,model=x$model)#new_INGARCH_EM(pb.dat,init=x$estimate_vector,model=x$model)  
        pb.cor[pb] = cor(pb.dat)[1,2]
        pb.mat[pb,] = pb.em$estimate_vector
        pb=pb+1
      },
      error = function(e){ 
      }
    )
    
    #pb.mean[pb,] = colMeans(pb.dat)
    #pb.var[pb,] = apply(pb.dat,2,var)
    #pb.kurtosis[pb,] = diag(kurtosis(pb.dat))
    #pb.skew[pb,] = apply(pb.dat,2,skewness)
    #pb.skew.nb[pb,] = pb.skew[pb,] - (2*pb.var[pb,]-pb.mean[pb,])/pb.mean[pb,]/sqrt(pb.var[pb,])
    #pb.quantile1[pb,] = c(summary(pb.dat[,1]))
    #pb.quantile2[pb,] = c(summary(pb.dat[,2]))
    if(pb%%10==0) {print(pb)}
  }
  ci.mat = apply(pb.mat,2,quantile,c(0.025,0.975))
  pb =list('95%_CI'= ci.mat, average = colMeans(pb.mat),pb.est = pb.mat,pb.cor=pb.cor)
  class(pb) = "ingarch_pb"
  return(pb)
}
print.ingarch_pb = function(x){print("use 'names'.")}







####Marginal PIT histograms
####Marginal PIT histograms


##Parameterization mu=R(phi,alpha)
pit_mSichel_old = function(dat, ingarch_em, series){ #Only bivariate. #ingarch_em = result of em.algorithm.
  #length(dat1)=n-max.b, length(lam)=n-max.b
  #dat1 = series of interest, dat2 = series to be summed out.
  max.ab = max(ingarch_em$model[["past_obs"]])
  dat1 = dat[-(1:max.ab),series]
  dat2 = dat[-(1:max.ab),-series]
  phi = ingarch_em$estimate$phi
  alpha = ingarch_em$estimate$alpha
  lambda = fitted(ingarch_em)
  
  
  #P_t (y_t)
  Py_t = rep(0,length(dat1));  Py_past = rep(0,length(dat1))
  for(i in 1:length(dat1)){
    tmp = 0
    for(j in 0:dat1[i]){
      grd=cbind(0:600,0:600)
      grd[,series] = j
      tmp.add = sum(apply( grd,1,function(x){ exp(old_mSichel_pairs(x,phi2=phi,alpha2=alpha,lambda2=lambda[i,])) }) )
      tmp = tmp + tmp.add
    }
    Py_t[i] = tmp
    if(dat1[i]>0){ Py_past[i] = tmp- tmp.add} 
    print(i)
  }
  return(list(Py_t=Py_t,Py_past=Py_past))
}


pit_u = function(u,pit_ingarch){  #pit_ingarch = list(Py_t,Py_past)
  Py_t = pit_ingarch[["Py_t"]]
  Py_past = pit_ingarch[["Py_past"]]
  F_t = rep(0,length(Py_t))
  F_t[Py_past>=u] = 0
  F_t[which(Py_t<=u)] = 1
  F_t[Py_t>u & Py_past<u] = ( (u-Py_past)/(Py_t-Py_past) )[Py_t>u & Py_past<u]
  return(mean(F_t))
}
pit_mar = Vectorize(pit_u ,vectorize.args = "u")





##Univariate PIT histogram.
##Univariate PIT histogram.
##Univariate PIT histogram.
Sichel_u = function(dat,ingarch_em,series){ #length(dat1)=n-max.b, length(lam)=n-max.b
  
  #P_t (y_t)
  max.ab = max(ingarch_em$model[["past_obs"]])
  dat1 = dat[-(1:max.ab),series]
  lam1 = fitted(ingarch_em)[,series]
  phi2 = ingarch_em$estimate$phi
  alpha2 =ingarch_em$estimate$alpha
  
  Py_t = apply(cbind(dat1,lam1),1,function(x){pSICHEL(q=x[1],mu=x[2],sigma=1/phi2,nu=alpha2)})
  Py_past = Py_t
  Py_past[dat1==0] = 0
  Py_past[dat1>0] = apply(cbind(dat1[dat1>0],lam1[dat1>0]),1,function(x){pSICHEL(q=x[1]-1,mu=x[2],sigma=1/phi2,nu=alpha2)} )
  
  return(list(Py_t=Py_t,Py_past=Py_past))
}



old_Sichel_u = function(dat,ingarch_em,series){ #length(dat1)=n-max.b, length(lam)=n-max.b
  
  #P_t (y_t)
  max.ab = max(ingarch_em$model[["past_obs"]])
  dat1 = dat[-(1:max.ab),series]
  lam1 = fitted(ingarch_em)[,series]
  phi2 = ingarch_em$estimate$phi
  alpha2 =ingarch_em$estimate$alpha
  
  Py_t = apply(cbind(dat1,lam1),1,function(x){pSICHEL(q=x[1],mu=x[2]*besselRatio(phi2,alpha2,1),sigma=1/phi2,nu=alpha2)})
  Py_past = Py_t
  Py_past[dat1==0] = 0
  Py_past[dat1>0] = apply(cbind(dat1[dat1>0],lam1[dat1>0]),1,function(x){pSICHEL(q=x[1]-1,mu=x[2]*besselRatio(phi2,alpha2,1),sigma=1/phi2,nu=alpha2)} )
  
  return(list(Py_t=Py_t,Py_past=Py_past))
}


Poi_u = function(dat,ingarch_em,series){ #length(dat1)=n-max.b, length(lam)=n-max.b
  
  max.ab = max(ingarch_em$model[["past_obs"]])
  dat1 = dat[-(1:max.ab),series]
  lam1 = fitted(ingarch_em)[,series]
  
  #P_t (y_t)
  Py_t = apply(cbind(dat1,lam1),1,function(x){ppois(q=x[1],lambda=x[2])})
  Py_past = Py_t
  Py_past[dat1==0] = 0
  Py_past[dat1>0] = apply(cbind(dat1[dat1>0],lam1[dat1>0]),1,function(x){ppois(q=x[1]-1,lambda=x[2])} )
  
  return(list(Py_t=Py_t,Py_past=Py_past))
}










####Functions for Poisson model (for comparison purpose)
neg.q.log.lik = function(dat,pars,model){   
  #cond_z = a vector of E(Z_t|y_t).
  
  max.ab = max(model[["past_obs"]])
  lambdam = mean_process(dat,pars,model)
  
  result = -sum(rowSums(lambdam)) + sum(rowSums(log(lambdam)*dat[-(1:max.ab),,drop=F]))
  return(-result)
}











###Poisson model (Fokianos, 2020.)
new_INGARCH_Poi = function(dat,init=NULL,model,method="BFGS"){
  
  init_val = init
  if( is.null(init) == T ){ init_val = new_INGARCH_init(dat,model,distr = "Poisson")$estimate}
  
  dat_n = dim(dat)[1] #the number of observation
  dat_p = dim(dat)[2] #the dimension of multivariate time series; 
  
  bb = length(model[["past_obs"]])
  aa = length(model[["past_mean"]])
  max.ab = max(model[["past_obs"]])
  tmp_dat = dat[-c(1:max.ab),]
  
  num_par = dat_p + (bb*dat_p^2 ) + (aa*dat_p^2 )
  lambda = mean_process(dat,init_val,model)
  
  new.est = 
    optim(init_val, function(pars){
      neg.q.log.lik(dat,pars,model)},method=method)
  new.est.par = new.est$par  
  lambda = mean_process(dat,new.est.par,model)
  if( is.null(model[["past_mean"]])==F ){ #GARCH
    results = list(estimate_vector=new.est,par, log_QL = -new.est$value, fitted=lambda, residuals = tmp_dat-lambda, #QIC = ,
                   estimate = list(d=new.est.par[1:dat_p],B=matrix(new.est.par[(1:(bb*dat_p^2))+dat_p],nrow=dat_p),
                                   A=matrix(new.est.par[(1:(aa*dat_p^2))+dat_p+bb*dat_p^2],nrow=dat_p)),sample_size=dat_n, num_par = num_par ,model=model)
  }else{ #ARCH
    results = list(estimate_vector=new.est,par, log_QL = -new.est$value, fitted=lambda, residuals = tmp_dat-lambda, #QIC = ,
                   estimate = list(d=new.est.par[1:dat_p],B=matrix(new.est.par[(1:(bb*dat_p^2))+dat_p],nrow=dat_p),A=NULL),
                   sample_size=dat_n, num_par = num_par, model=model)  
  }    
  class(results) = "ingarch"
  return(results)
}
print.ingarch = function(x){print(x[c("estimate","log_QL")])}




