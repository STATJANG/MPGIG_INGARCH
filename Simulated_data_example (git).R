###Simulation Setting
load("MPGIG_INGARCH_FN (git).R")
samplesize=1000 ; mc_size = 500 #Repetition size. 

#True parameters
true_phi = 0.5; true_alpha= -1.5
true_d = c(1,1);
dat_p = length(true_d) 
true_A = matrix(c(0.3,0,0,0.25),ncol=dat_p)
true_B = matrix(c(0.4,0,0,0.3),ncol=dat_p)
nustart=rep(0,dat_p)

#Stationarity and ergodicity conditions
sqrt(max(abs(eigen(t(true_B)%*%true_B)[[1]])))+sqrt(max(abs(eigen(t(true_A)%*%true_A)[[1]])))
max(colSums(abs(true_A)))


#indexing banded diagonal elements if needed.
band_index = abs(col(true_A)-row(true_A))<=1
diag_index = (abs(col(true_A)-row(true_A))==0)

num_par = 2 + dat_p + 2*dat_p^2 #number of pars
stored_est = matrix(0,nrow=mc_size,ncol=num_par)

     set.seed(1)
     for(simul_mc in 1:mc_size){
          #Data generation
          datt = sim_pgig_ingarch(samplesize,phi5 = true_phi,alpha5 = true_alpha,
                                  d5=true_d, A5 = true_A, B5 = true_B, burns = 500) 
          #Parameter estimation
          tmp_val = old_INGARCH_EM(datt,model=list(past_mean=1,past_obs=1),trace=T)
          
          stored_est[simul_mc,] = tmp_val$estimate_vector
          cat("simul_mc=",simul_mc)
          }
