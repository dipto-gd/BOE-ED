rm(list=ls())

#install.packages("StanHeaders", lib = "/home/faird/ghosh189/R/x86_64-pc-linux-gnu-library/4.3/")
#install.packages("RcppParallel", lib = "/home/faird/ghosh189/R/x86_64-pc-linux-gnu-library/4.3/")


## MSI locations
#library(parallel)
.libPaths("/home/faird/ghosh189/R/x86_64-pc-linux-gnu-library/4.3/")
library("StanHeaders")
library("rstan")

library(MASS)
library(dplyr)
# if (!require(mvtnorm)) {
#   install.packages("mvtnorm")
# }
library(mvtnorm,lib.loc="/home/faird/ghosh189/R/x86_64-pc-linux-gnu-library/4.0")
setwd("/panfs/jay/groups/4/ventz001/ghosh189/Opt_adap_design/New_Codes")
source("imp_functions.R")


## local locations
# source("/Users/diptogd/Documents/Academic /Research/NLP_Codes/imp_functions.R")


#setwd("/panfs/jay/groups/4/ventz001/ghosh189/Opt_adap_design/New_Codes")

#setwd("/panfs/jay/groups/4/ventz001/ghosh189/Opt_adap_design/New_Codes/results")
#setwd("C:/Users/User/Desktop/Academic/UMN/Fall 2023/Steffen-RA/Simulations/Local Codes/local_results")

no_scenarios = 17


N_trials = 1000
no_jobs = 50

n_stage1_trials = round(N_trials/no_jobs)

Job_mat = expand.grid(Job=1:no_jobs,Sc=1:(no_scenarios))

# job runs from no_jobs *no_scenarios
a = as.integer(commandArgs()[3])

print(a)


time1       = Sys.time()
rho1        = 1
rho2        = 0
lambda      = 5
num_cores   = 10

indices     = c(1:4)
n.subsets   = length(indices)
subsets     = lapply(1:(2^n.subsets - 1), function(i) indices[which(bitwAnd(i, 2^(0:(n.subsets - 1))) != 0)])
subsets.string =  sapply(subsets, function(x) paste(x, collapse = ","))


Sc = Job_mat[a,"Sc"]
Batch = Job_mat[a,"Job"]
set.seed(a) 
time1=Sys.time()

prop     = scenarios_prop_phi_nlp(Sc)$prop
prop_ED  = scenarios_prop_phi_nlp(Sc)$prop_ED
b.phi      = scenarios_prop_phi_nlp(Sc)$b.phi
b.phi_ED   =  scenarios_prop_phi_nlp(Sc)$b.phi_ED
intercept = scenarios_prop_phi_nlp(Sc)$intercept

bx345 = c(0.5,0.5,0.5) 
b.phi = c(b.phi,bx345)
b.phi_ED = c(b.phi_ED,bx345)

# no of pretrt covariates

prob_T_Rct = 2/3
prob_T_Ed  = 1/20

phi        = b.phi[3:6]
phi_ED     = b.phi_ED[3:6]
size.subsets           = lapply(subsets,function(indices.vec) sum(prop[indices.vec])^rho1)

n_precov = length(b.phi_ED)-length(phi)

weighted.phi.H0.true         = sapply(subsets,function(indices.vec)(prop[indices.vec]) %*% phi[indices.vec] /sum(prop[indices.vec])) 

N1        = 100
N2        = 100
N_ed      = 5000
v.Y       = 1

# scalar for local/ nonlocal prior
s         = 0.05
s_intercept = 0.5
s_lp        =1
intercept_0 = intercept
constant    = rep(-0.1,length(subsets))   
futil_thresh = 0.2 ###################### similar to other futlity stopping criteria

prior.Ivar = diag(1/100, length(b.phi))
prior.mu   = rep(0,length(b.phi))



#' ============ @data-stage-1
time1=Sys.time()

# Generate n_stage1_trials 
# n_stage1_trials = 10000

#n_stage1_trials = 10

data_stage1_list=replicate(n = n_stage1_trials , simplify=F, expr = {
  
  #'===========@internal_data
  indices= c(1,2,3,4)
  
  R.i_RCT       = data_gen2(indices,N1,prop,prob_T_Rct)
  
  Y.i_RCT    =  R.i_RCT%*%b.phi + rnorm(N1)
  
  #prop of control is given by
  # mean(apply(R.i, 1, function(row) all(row[3:6] == 0)))
  
  #'===========@external_data
  indices= c(1,2,3,4)
  
  R.i_ED       = data_gen2(indices,N_ed,prop_ED,prob_T_Ed)
  
  Y.i_ED    = intercept+ R.i_ED%*%b.phi_ED + rnorm(N_ed)
  
  R.i       = rbind(R.i_RCT,R.i_ED)
  Y.i       = rbind(Y.i_RCT,Y.i_ED)
  list(Y.i.H0=Y.i, R.i=R.i)
})



#' ============ @loop_over_stage1_datasets
setwd("/panfs/jay/groups/4/ventz001/ghosh189/Opt_adap_design/Rstan_codes")

# #When ED intercept has local prior, gamma also has local prior
# c_code_localpriorED_alt= stan_model(file='lin_mod_ED_intercept_lp.stan',verbose = T, save_dso = T, auto_write = T)
# 
# #When ED may be different, using nonLocal prior
# c_code_nlppriorED_alt= stan_model(file='lin_mod_ED_intercept_nlp.stan',verbose = T, save_dso = T, auto_write = T)

#setwd("/panfs/jay/groups/4/ventz001/ghosh189/Opt_adap_design/Rstan_codes")

# When all gammas are known to be 0 (Null)
c_code_noprior_null= stan_model(file='lin_mod_null.stan',verbose = T, save_dso = T, auto_write = T)

# When all 4 gammas have a local prior and  there is no intercept
c_code_localprior_for_4interactions= stan_model(file='lin_mod_localprior_for_4interactions.stan',verbose = T, save_dso = T, auto_write = T)

# When ED intercept has a local prior and no gammas
c_code_lppriorED_alt_no_gamma= stan_model(file='lin_mod_ED_intercept_lp_gamma_null.stan',verbose = T, save_dso = T, auto_write = T)

# When ED intercept has a non local prior and no gammas
c_code_nlppriorED_alt_no_gamma= stan_model(file='lin_mod_ED_intercept_nlp_gamma_null.stan',verbose = T, save_dso = T, auto_write = T)


### New NLP with bias_beta codes

#When ED intercept, betabias have local prior, gammabias has local prior, gamma also has local prior
c_code_localpriorED_alt= stan_model(file='lin_mod_ED_intercept_biasbeta_lp.stan',verbose = T, save_dso = T, auto_write = T)

#When ED may be different, using nonLocal prior for intercept, gammabias and betabiashas local prior,
c_code_nlppriorED_alt= stan_model(file='lin_mod_ED_biasbeta_intercept_nlp.stan',verbose = T, save_dso = T, auto_write = T)

# # When ED intercept and beta bias have local priors and no gammas
# c_code_lppriorED_alt_no_gamma = stan_model(file='ED_betabias_lp_intercept_lp_no_gamma.stan',verbose = T, save_dso = T, auto_write = T)
# 
# # When ED intercept has nlp and beta bias has local priors and no gammas
# c_code_nlppriorED_alt_no_gamma = stan_model(file='ED_betabias_lp_intercept_nlp_no_gamma.stan',verbose = T, save_dso = T, auto_write = T)


M = expand.grid(0:1, 0:1, 0:1, 0:1)
D = expand.grid(0:1)

#' following code lines give optimal decision at the end of stage 1

#time1=Sys.time()

log_liks_stage1 = function(R.i, Y.i, label.i) lapply(1:nrow(D), function(j){
  trt_int   = R.i[,c(3:6)]
  if(j==1){
    
    
    logprob_data_given_models_lp_nlp =lapply(1:nrow(M), function(model_num){
      m = M[model_num,]
      i.m = which(m==1)
      G = sum(m==1)
      
      trt_int = R.i[, c(3:6)[i.m]] # selects only those columns that cause interactions with trt among columns 3:6 under model m
      
      
      
      
      if (model_num==1){
        #'====================@null-no-gamma_no_ED
        
        # samples_nop = sampling(object = c_code_noprior_null ,data = data)
        # samples_MCMC_nop=extract(samples_nop)
        
        data = list(y=as.vector(Y.i),
                    x= R.i[,c(1,2,7:9)],
                    K=5,
                    G=G,
                    s=1,
                    s_lp=1,
                    N=nrow(R.i))
        opt_result_nop = optimizing(c_code_noprior_null,hessian=T,  data = data,importance_resampling=TRUE,draws=4000)
        
        importance_samples = opt_result_nop$theta_tilde 
        log_importance_weights_unscaled = opt_result_nop$log_p-opt_result_nop$log_g
        ### Importance weights log_p-log_g are just ratios of densities not scaled by normalization constant 
        log_importance_weights = log_importance_weights_unscaled
        #range(log_importance_weights)
        
        posterior_mode_without_bias_nop = opt_result_nop$par
        
        hessian_matrix_without_bias_nop = opt_result_nop$hessian
        
        posterior_mode =posterior_mode_without_bias_nop
        hessian_matrix= hessian_matrix_without_bias_nop
        
        beta_hat   = posterior_mode_without_bias_nop[grep("^beta",names(posterior_mode_without_bias_nop))]
        gamma_hat  = rep(0,4)
        gamma_var  = matrix(0,nrow=4,ncol=4)
        mu_hat     = posterior_mode_without_bias_nop["mu"]
        sigma_hat  = posterior_mode_without_bias_nop["sigma"]
        likelihood_data = sum(dnorm(as.vector(Y.i) , mean = R.i[,c(1,2,7:9)]%*% beta_hat , sd=sigma_hat,log=T)) +
          sum(dnorm(beta_hat, mean= mu_hat, 1,log=T)) +dnorm(mu_hat,0,1,log=T) + 
          dnorm(sigma_hat,0,5,log=T)
        
        logprob_data_given_model =  -0.5* determinant(-1*hessian_matrix_without_bias_nop,logarithm = T)$modulus +  0.5*length(posterior_mode_without_bias_nop)*log(2*pi) +likelihood_data
        
      }else{
        #'====================@local-prior
        # samples = sampling(object = c_code_localprior_for_4interactions ,data = data)
        # samples_MCMC=extract(samples)
        data = list(y=as.vector(Y.i),
                    x= R.i[,c(1,2,7:9)],
                    trt_int = as.matrix(trt_int),
                    K=5,
                    G=G,
                    s=s,
                    s_lp=1,
                    N=nrow(R.i))
        
        opt_result = optimizing(c_code_localprior_for_4interactions,hessian=T,  data = data,importance_resampling=TRUE,draws=4000)
        
        importance_samples = opt_result$theta_tilde
        log_importance_weights_unscaled = opt_result$log_p-opt_result$log_g
        ### Importance weights log_p-log_g are just ratios of densities not scaled by normalization constant
        log_importance_weights = log_importance_weights_unscaled
        #range(log_importance_weights)
        
        posterior_mode_localprior = opt_result$par
        
        hessian_matrix_localprior = opt_result$hessian
        
        posterior_mode =posterior_mode_localprior
        hessian_matrix= hessian_matrix_localprior
        
        beta_hat   = posterior_mode_localprior[grep("^beta",names(posterior_mode_localprior))]
        gamma_hat  =  posterior_mode_localprior[grep("^gamma",names(posterior_mode_localprior))]
        mu_hat     = posterior_mode_localprior["mu"]
        sigma_hat  = posterior_mode_localprior["sigma"]
        likelihood_data = sum(dnorm(as.vector(Y.i) , mean = R.i[,c(1,2,7:9)]%*% beta_hat + as.matrix(trt_int)%*% gamma_hat, sd=sigma_hat,log=T)) +
          sum(dnorm(beta_hat, mean= mu_hat, 1,log=T)) +dnorm(mu_hat,0,1,log=T) +
          sum(dnorm(gamma_hat, mean= 0, 1,log=T))+ dnorm(sigma_hat,0,5,log=T)
        
        
        logprob_data_given_model =  -0.5* determinant(-1*hessian_matrix_localprior,logarithm = T)$modulus +  0.5*length(posterior_mode_localprior)*log(2*pi) +likelihood_data
        
        
      }
      return(list(logprob_data_given_model=logprob_data_given_model,
                  posterior_mode=posterior_mode,
                  hessian_matrix=hessian_matrix,
                  importance_samples=importance_samples,
                  log_importance_weights=log_importance_weights))
    })
    
    
  }else if(j==2){
    
    logprob_data_given_models_lp_nlp =lapply(1:nrow(M), function(model_num){
      m = M[model_num,]
      i.m = which(m==1)
      G = sum(m==1)
      
      trt_int = R.i[, c(3:6)[i.m]] # selects only those columns that cause interactions with trt among columns 3:6 under model m
      
      
      
      
      if (model_num==1){
        #'====================@null
        
        # samples_nop = sampling(object = c_code_noprior_null ,data = data)
        # samples_MCMC_nop=extract(samples_nop)
        
        data = list(y=as.vector(Y.i),
                    x= R.i[,c(1,2,7:9)],
                    label=label.i,
                    K=5,
                    G=G,
                    s=s,
                    s_intercept=s_intercept,
                    intercept_0= intercept_0,
                    s_lp=s_lp,
                    N=nrow(R.i))
        #'====================@local-prior
        
        opt_result = optimizing(c_code_lppriorED_alt_no_gamma,hessian=T,  data = data)
        
        posterior_mode_localprior = opt_result$par
        
        #'====================@nonlocal-prior
        
        opt_result_nop = optimizing(c_code_nlppriorED_alt_no_gamma,hessian=T, init=as.list(posterior_mode_localprior), data = data,importance_resampling=TRUE,draws=4000)
        
        importance_samples = opt_result_nop$theta_tilde 
        log_importance_weights_unscaled = opt_result_nop$log_p-opt_result_nop$log_g
        
        ### Importance weights log_p-log_g are just ratios of densities not scaled by normalization constant 
        log_importance_weights = log_importance_weights_unscaled
        
        posterior_mode_without_bias_nop = opt_result_nop$par
        
        hessian_matrix_without_bias_nop = opt_result_nop$hessian
        
        posterior_mode =posterior_mode_without_bias_nop
        hessian_matrix= hessian_matrix_without_bias_nop
        
        beta_hat   = posterior_mode_without_bias_nop[grep("^beta",names(posterior_mode_without_bias_nop))]
        intercept_hat = posterior_mode_without_bias_nop["intercept"]
        gamma_hat  = rep(0,4)
        gamma_var  = matrix(0,nrow=4,ncol=4)
        mu_hat     = posterior_mode_without_bias_nop["mu"]
        sigma_hat  = posterior_mode_without_bias_nop["sigma"]
        likelihood_data = sum(dnorm(as.vector(Y.i) , mean = R.i[,c(1,2,7:9)]%*% beta_hat+intercept_hat*label.i , sd=sigma_hat,log=T)) +
          sum(dnorm(beta_hat, mean= mu_hat, 1,log=T)) +dnorm(mu_hat,0,1,log=T) + 
          sum(dnorm(intercept_hat,mean=0,s_intercept,log=T))+log(intercept_hat^2) - length(intercept_hat)*log(s_intercept^2) +
          dnorm(sigma_hat,0,5,log=T)
        
        
        logprob_data_given_model =  -0.5* determinant(-1*hessian_matrix_without_bias_nop,logarithm = T)$modulus +  0.5*length(posterior_mode_without_bias_nop)*log(2*pi) +likelihood_data
        
        
        
      }else{
        data = list(y=as.vector(Y.i),
                    x= R.i[,c(1,2,7:9)],
                    trt_int = as.matrix(trt_int),
                    label=label.i,
                    K=5,
                    G=G,
                    s=1,
                    s_intercept=s_intercept,
                    intercept_0= intercept,
                    s_lp=s_lp,
                    N=nrow(R.i))
        
        opt_result = optimizing(c_code_localpriorED_alt,hessian=T,  data = data,importance_resampling=TRUE,draws=4000)
        
        posterior_mode_localprior = opt_result$par
        
        
        #'====================@nonlocal-prior
        
        # samples_nlp = sampling(object = c_code_nonlocalprior_for_4interactions ,data = data)
        # samples_MCMC_nlp=extract(samples_nlp)
        
        
        opt_result_nlp = optimizing(c_code_nlppriorED_alt,hessian=T, init=as.list(posterior_mode_localprior),data = data,importance_resampling=TRUE,draws=4000)
        
        
        importance_samples = opt_result_nlp$theta_tilde 
        log_importance_weights_unscaled = opt_result_nlp$log_p-opt_result_nlp$log_g
        ### Importance weights log_p-log_g are just ratios of densities not scaled by normalization constant 
        log_importance_weights = log_importance_weights_unscaled
        
        posterior_mode_nlp = opt_result_nlp$par
        
        hessian_matrix_nlp = opt_result_nlp$hessian
        
        posterior_mode =posterior_mode_nlp
        hessian_matrix= hessian_matrix_nlp
        
        beta_hat   = posterior_mode_nlp[grep("^beta",names(posterior_mode_nlp))]
        intercept_hat = posterior_mode_nlp["intercept"]
        gamma_hat  =  posterior_mode_nlp[grep("^gamma",names(posterior_mode_nlp))]
        mu_hat     = posterior_mode_nlp["mu"]
        sigma_hat  = posterior_mode_nlp["sigma"]
        
        likelihood_data = sum(dnorm(as.vector(Y.i) , mean = R.i[,c(1,2,7:9)]%*% beta_hat + as.matrix(trt_int)%*% gamma_hat+intercept_hat*label.i, sd=sigma_hat,log=T)) +
          sum(dnorm(beta_hat, mean= mu_hat, 1,log=T)) +dnorm(mu_hat,0,1,log=T) +
          sum(dnorm(intercept_hat,mean=0,s_intercept,log=T))+ log(intercept_hat^2) - length(intercept_hat)*log(s_intercept^2)+
          sum(dnorm(gamma_hat, mean= 0, 1,log=T))+ log(2)+dnorm(sigma_hat, mean=0,5, log=T)
        
        logprob_data_given_model =  -0.5* determinant(-1*hessian_matrix_nlp,logarithm = T)$modulus +  0.5*length(posterior_mode_nlp)*log(2*pi) +likelihood_data
        
        
      }
      return(list(logprob_data_given_model=logprob_data_given_model,
                  posterior_mode=posterior_mode,
                  hessian_matrix=hessian_matrix,
                  importance_samples=importance_samples,
                  log_importance_weights=log_importance_weights))
    })
    
  }
})

posterior_gamma_corrected = function(logprob_data_given_models_lp_nlp_total_unlisted) lapply(1:length(logprob_data_given_models_lp_nlp_total_unlisted), function(iter){
  label = iter
  posterior_mode_vec =  logprob_data_given_models_lp_nlp_total_unlisted[[label]]$posterior_mode
  
  posterior_mode = matrix(0,nrow=1,ncol= length(b.phi)+3)
  posterior_mode[1:n_precov] =posterior_mode_vec[1:n_precov]
  model_num = ifelse(label>16,label-16,label)
  
  param_names =names(posterior_mode_vec)
  id_cols =grep("gamma",param_names)
  if(length(id_cols)>0){
    
    posterior_mode[which(M[model_num,]!=0)+n_precov]= posterior_mode_vec[id_cols]
  }
  id_cols =grep("intercept",param_names)
  if(length(id_cols)>0){
    
    posterior_mode[length(posterior_mode)]  = posterior_mode_vec["intercept"]
  }
  posterior_mode[length(posterior_mode)-2]=posterior_mode_vec["sigma"]
  posterior_mode[length(posterior_mode)-1]=posterior_mode_vec["mu"]
  
  colnames(posterior_mode) =c(param_names[1:n_precov], paste0("gamma[",1:length(phi),"]"),"sigma","mu","intercept")
  
  return(posterior_mode)
})
time1 = Sys.time()
data_stage_12_full_list = sapply(1:n_stage1_trials, function(i){
  
  D1=data_stage1_list[[i]]
  
  R.i       = D1$R.i
  Y.i       = D1$Y.i.H0
  label.i=c(rep(0,N1),rep(1,N_ed))
  
  #G         = 4
  logprob_data_given_models_lp_nlp_total = log_liks_stage1(R.i, Y.i, label.i)
  
  
  logprob_data_given_models_lp_nlp_total_unlisted = do.call(c, logprob_data_given_models_lp_nlp_total)
  
  
  
  ## This is now a list of lists.
  logliks= sapply(logprob_data_given_models_lp_nlp_total_unlisted, function(V)V$logprob_data_given_model)
  which.max(logliks)
  model_wts = sapply(logliks, function(ll){1/(sum(exp(logliks-ll)))})
  plot(1:32, model_wts)
  
  ## Returns reformatted posterior modes to incorporate 0 gammas whenever they are 0
  posterior_mode_list = posterior_gamma_corrected(logprob_data_given_models_lp_nlp_total_unlisted)
  
  ## calculating u1 
  unweighted_u1 = sapply(1:length(logprob_data_given_models_lp_nlp_total_unlisted), function(iter){ 
    
    posterior_mode = posterior_mode_list[[iter]]
    gamma = posterior_mode[,n_precov+1:length(phi)]
    util_subsets_func(gamma ,subsets,prop,1,0)
    
    
  })
  
  weighted_u1 = unweighted_u1 %*% model_wts
  
  
  ## Draw model indicators M_i from these posterior model probs n_iter =10000 times
  n_iter  = 1000
  samples = sample.int(length(model_wts), size = n_iter,replace = T, prob = model_wts)
  
  ## Draw theta_i from the model posteriors, get the weights of these theta_i's
  sampled_models = table(samples)
  Theta_i_pred_list =lapply(1:length(sampled_models), function(iter){
    label = as.numeric(names(sampled_models)[iter])
    freq = sampled_models[iter]
    importance_samples =  logprob_data_given_models_lp_nlp_total_unlisted[[label]]$importance_samples
    
    draws = matrix(0,nrow=nrow(importance_samples),ncol= length(b.phi)+3)
    draws[,1:n_precov] =importance_samples[,1:n_precov]
    model_num = ifelse(label>16,label-16,label)
    
    param_names =colnames(importance_samples)
    id_cols =grep("gamma",param_names)
    if(length(id_cols)>0){
      
      draws[,which(M[model_num,]!=0)+n_precov]= importance_samples[,id_cols]
    }
    id_cols =grep("intercept",param_names)
    if(length(id_cols)>0){
      
      draws[,ncol(draws)]  = importance_samples[,"intercept"]
    }
    draws[,ncol(draws)-2]=importance_samples[,"sigma"]
    draws[,ncol(draws)-1]=importance_samples[,"mu"]
    
    colnames(draws) =c(colnames(importance_samples)[1:n_precov], paste0("gamma[",1:length(phi),"]"),"sigma","mu","intercept")
    
    
    log_importance_weights = logprob_data_given_models_lp_nlp_total_unlisted[[label]]$log_importance_weights
    log_importance_weights = log_importance_weights-min(log_importance_weights)
    importance_weights = exp(log_importance_weights)
    posterior_samples_indices = sample(x =nrow(draws), size= freq, rep=T, prob =importance_weights)
    
    
    if(freq==1){
      posterior_samples=matrix(c(label,draws[posterior_samples_indices,]),nrow=1)
    }else{
      posterior_samples=cbind(label,draws[posterior_samples_indices,])
    }
    #sampled_log_importance_weights = log_importance_weights[posterior_samples_indices]
    #hessian_matrix     = logprob_data_given_models_lp_nlp_total_unlisted[[label]]$hessian_matrix
    return(posterior_samples)
    #sampled_log_importance_weights=sampled_log_importance_weights,
    #hessian_matrix=hessian_matrix 
    
    #rmvnorm(1,mean= posterior_mode, hessian_matrix)
  })
  
  
  ## make an array of posterior theta_i by splitting rows when multiple thetas
  M_i_Theta_i_pred_list = NULL
  for(z in Theta_i_pred_list) M_i_Theta_i_pred_list = rbind(M_i_Theta_i_pred_list,z)
  
  ## append the columns of 16 CATES for each theta in M_i_Theta_i_pred_list
  
  M_i_Theta_i_phi_i_pred_list = matrix(NA,nrow=nrow(M_i_Theta_i_pred_list),ncol=ncol(M_i_Theta_i_pred_list)+length(subsets))
  for(z in 1:nrow(M_i_Theta_i_pred_list)){
    phi  = M_i_Theta_i_pred_list[z, (2+n_precov):(5+n_precov)]
    M_i_Theta_i_phi_i_pred_list[z,] = c(M_i_Theta_i_pred_list[z,],sapply(subsets,function(indices.vec)(prop[indices.vec]) %*% phi[indices.vec] /sum(prop[indices.vec])) )
  } 
  colnames(M_i_Theta_i_phi_i_pred_list) =c(colnames(M_i_Theta_i_pred_list),paste0("phi_",1:15))
  prob_futil = max(colMeans( M_i_Theta_i_phi_i_pred_list[,13+1:length(subsets)]>0))
  u0 = ifelse(prob_futil>futil_thresh,-1000,1000)
  ## Draw X_i from original true population
  
  X_i_2_list   = lapply(1:n_iter,function(model_num){
    X_i_2 = lapply(1:4,function(i)data_gen2(i,N2,prop,prob_T_Rct))
    X_i_2_mat = do.call(rbind, X_i_2)
  })
  ####################################################################
  # Make sure to draw the 4 unique subgroups exactly 100 time each
  ## Add intercept_hat in Theta_i matrix for the models
  
  ## down below make intercept =0 if absent, and get rid of ifelse
  
  ## after the thetas, get the cate for each of 1000 rows
  ####################################################################
  
  ## Given X_i, theta_i, draw Y_i_2
  
  Y_i_2 = sapply(1:n_iter, function(iter){
    
    
    posterior_theta_vec = M_i_Theta_i_phi_i_pred_list[iter,]
    model_num           = posterior_theta_vec[[1]]
    gamma_model         = ifelse(model_num>=17,model_num-nrow(M),model_num) 
    param_names         = names(posterior_theta_vec)
    beta_hat            = posterior_theta_vec[grep("^beta\\[", param_names)]
    gamma_hat           = posterior_theta_vec[grep("^gamma\\[", param_names)]
    sigma_hat           = posterior_theta_vec["sigma"]
    intercept_hat       = posterior_theta_vec["intercept"]
    b.phi_pos = c(beta_hat[1:2],gamma_hat,beta_hat[3:5])
    Y.i.2_RCT    =  X_i_2_list[[iter]]%*%b.phi_pos + rnorm(n=nrow(X_i_2_list[[iter]]),0,sigma_hat)
    
    Y.i.2_RCT
  })
  
  
  utili_2_s = lapply(1:n_iter, function(iter_i) {
    
    ## iter_i is the outer loop running from 1 to n_iter.
    ## Finds the numerator term L(theta_i|D_2)
    
    #########################################################
    
    ## loop over subsets...
    ## sapply(subsets,
    
    sapply (1:length(subsets), function(subgp){
      
      indices.vec  = subsets[[subgp]]
      R.i.2         = X_i_2_list[[iter_i]]
      Y.i.2         = Y_i_2[,iter_i]
      
      ##... get subgp X matrix by filtering B,M 
      
      label_rep = if(length(indices.vec)==1){rep(indices.vec,N2)} else {sample(indices.vec,N2, replace=T,prop[indices.vec]/sum(prop[indices.vec]))}
      label_freq = table(label_rep)
      
      X_subset_rows =  as.vector(sapply(indices.vec, function(temp) ((temp-1)*100+1):(temp*100)))
      samples =unlist(sapply(indices.vec, function(temp){sample( ((temp-1)*100+1):(temp*100), label_freq[paste0(temp)],replace=FALSE)}))
      R_i_2_filt = R.i.2[samples,]
      
      ## get the corresponding Y's 
      Y_i_2_filt = Y.i.2[samples]
      
      ## get the likelihood weights for this particular D2|theta_i
      ## Inner loop runs also from 1 to n_iter
      log_likelihood_theta_1_c_given_D2 = sapply(1:n_iter, function(iter_j) {
        
        posterior_theta_vec = M_i_Theta_i_phi_i_pred_list[iter_j,]
        
        model_num           = posterior_theta_vec[1]
        param_names         = names(posterior_theta_vec)
        beta_hat            = posterior_theta_vec[grep("^beta\\[", param_names)]
        gamma_hat           = posterior_theta_vec[grep("^gamma\\[", param_names)]
        sigma_hat           = posterior_theta_vec[grep("^sigma", param_names)]
        
        sum(dnorm(as.vector(Y_i_2_filt), mean = R_i_2_filt[, c(1, 2, 7:9)] %*% beta_hat + R_i_2_filt[, c(3:6)] %*% gamma_hat, sd = sigma_hat, log = TRUE)) 
        
      })
      
      ## normalized EXPONENTIAL likelihood weights should be vector of 1000
      
      ##???????????????? check the next block of ??????
      normalized_likelihood_weights = exp(log_likelihood_theta_1_c_given_D2-max(log_likelihood_theta_1_c_given_D2))
      
      normalized_likelihood_weights= normalized_likelihood_weights/sum(normalized_likelihood_weights)
      
      #plot sorted normalized weights. check that the likeliood of the data from stage 2 is maximized at the parameter that generates it.
      #plot(1:1000,sort(normalized_likelihood_weights))
      
      ### may be use the normalized weights to get effective sample size 1/sum(normalized_wts^2) to reweigh utilities later
      
      #Sum (phi_j * normalized likelihood_weights)
      phi_subset = as.numeric(crossprod( M_i_Theta_i_phi_i_pred_list[,13+subgp],  normalized_likelihood_weights))
      phi_subset_tru = M_i_Theta_i_phi_i_pred_list[iter_i,13+subgp]
      
      
      ## D2*= I(e(phi_j|d2,d1)>-C?)
      d2_star = ifelse(phi_subset>-constant[subgp],1,0)
      ## cALCULATE U2 = U(D2*|first row theta i.e, theta[iter_i,])
      
      u2_star = ifelse(d2_star==0, ifelse(phi_subset_tru>0,0,-phi_subset_tru),ifelse(phi_subset_tru>0,phi_subset_tru+constant[subgp],constant[subgp]))
      
      ## return u2 
      u2_star
      ## sapply closes
      
    })
    ## so you have for each D2, 16 U2'S ie a 1000*16 matrix
    
    
    
  })
  
  
  ## col mean of u2 gives E(u2(d2*)|D1,d1=j) j=1,..16
  U2 = sapply(1:length(subsets),function(i)mean(sapply(utili_2_s,function(arr)arr[i])))
  
  total_U = U2+weighted_u1
  
  total_U = c(total_U, u0)
  d1_star = which.max(total_U)
  
  
  if(d1_star!=(length(subsets)+1)){
    # Simulate N2 stage 2 data points but only resticted to d1_star
    
    R_i_2  = data_gen2(subsets[[d1_star]],N2,prop,prob_T_Rct)
    
    
    ## Given X_i, draw Y_i_2 using true b.phi
    
    Y.i.2    = R_i_2 %*%b.phi + rnorm(N2)
    
    
    R.i.tot       = rbind(R.i,R_i_2)
    Y.i.tot       = rbind(Y.i,Y.i.2)
    label.i.tot  = c(label.i,rep(0,N2))
    
    logprob_data_given_models_lp_nlp_total = log_liks_stage1(R.i.tot, Y.i.tot, label.i.tot)
    
    logprob_data_given_models_lp_nlp_total_unlisted = do.call(c, logprob_data_given_models_lp_nlp_total)
    
    
    
    ## This is now a list of lists.
    logliks= sapply(logprob_data_given_models_lp_nlp_total_unlisted, function(V)V$logprob_data_given_model)
    which.max(logliks)
    model_wts = sapply(logliks, function(ll){1/(sum(exp(logliks-ll)))})
    
    ## Returns reformatted posterior modes to incorporate 0 gammas whenever they are 0
    posterior_mode_list = posterior_gamma_corrected(logprob_data_given_models_lp_nlp_total_unlisted)
    
    
    
    ## Draw model indicators M_i from these posterior model probs n_iter =10000 times
    n_iter  = 200
    samples = sample.int(length(model_wts), size = n_iter,replace = T, prob = model_wts)
    
    ## Draw theta_i from the model posteriors, get the weights of these theta_i's
    sampled_models = table(samples)
    Theta_i_pred_list =lapply(1:length(sampled_models), function(iter){
      label = as.numeric(names(sampled_models)[iter])
      freq = sampled_models[iter]
      importance_samples =  logprob_data_given_models_lp_nlp_total_unlisted[[label]]$importance_samples
      
      draws = matrix(0,nrow=nrow(importance_samples),ncol= length(b.phi)+3)
      draws[,1:n_precov] =importance_samples[,1:n_precov]
      model_num = ifelse(label>16,label-16,label)
      
      param_names =colnames(importance_samples)
      id_cols =grep("gamma",param_names)
      if(length(id_cols)>0){
        
        draws[,which(M[model_num,]!=0)+n_precov]= importance_samples[,id_cols]
      }
      id_cols =grep("intercept",param_names)
      if(length(id_cols)>0){
        
        draws[,ncol(draws)]  = importance_samples[,"intercept"]
      }
      draws[,ncol(draws)-2]=importance_samples[,"sigma"]
      draws[,ncol(draws)-1]=importance_samples[,"mu"]
      
      colnames(draws) =c(colnames(importance_samples)[1:n_precov], paste0("gamma[",1:length(phi),"]"),"sigma","mu","intercept")
      
      
      log_importance_weights = logprob_data_given_models_lp_nlp_total_unlisted[[label]]$log_importance_weights
      log_importance_weights = log_importance_weights-min(log_importance_weights)
      importance_weights = exp(log_importance_weights)
      posterior_samples_indices = sample(x =nrow(draws), size= freq, rep=T, prob =importance_weights)
      
      
      if(freq==1){
        posterior_samples=matrix(c(label,draws[posterior_samples_indices,]),nrow=1)
      }else{
        posterior_samples=cbind(label,draws[posterior_samples_indices,])
      }
      #sampled_log_importance_weights = log_importance_weights[posterior_samples_indices]
      #hessian_matrix     = logprob_data_given_models_lp_nlp_total_unlisted[[label]]$hessian_matrix
      return(posterior_samples)
      #sampled_log_importance_weights=sampled_log_importance_weights,
      #hessian_matrix=hessian_matrix 
      
      #rmvnorm(1,mean= posterior_mode, hessian_matrix)
    })
    
    
    ## make an array of posterior theta_i by splitting rows when multiple thetas
    M_i_Theta_i_pred_list = NULL
    for(z in Theta_i_pred_list) M_i_Theta_i_pred_list = rbind(M_i_Theta_i_pred_list,z)
    
    ## append the columns of 16 CATES for each theta in M_i_Theta_i_pred_list
    
    M_i_Theta_i_phi_i_pred_list = matrix(NA,nrow=nrow(M_i_Theta_i_pred_list),ncol=ncol(M_i_Theta_i_pred_list)+length(subsets))
    for(z in 1:nrow(M_i_Theta_i_pred_list)){
      phi  = M_i_Theta_i_pred_list[z, (2+n_precov):(5+n_precov)]
      M_i_Theta_i_phi_i_pred_list[z,] = c(M_i_Theta_i_pred_list[z,],sapply(subsets,function(indices.vec)(prop[indices.vec]) %*% phi[indices.vec] /sum(prop[indices.vec])) )
    } 
    colnames(M_i_Theta_i_phi_i_pred_list) =c(colnames(M_i_Theta_i_pred_list),paste0("phi_",1:15))
    
    phi_j_d1_star = mean(M_i_Theta_i_phi_i_pred_list[,paste0("phi_",d1_star)])
    ## D2*= I(e(phi_j|d2,d1)>-C?)
    d2_star = ifelse(phi_j_d1_star>-constant[d1_star],1,0)
  }else{d2_star =0 }
  c(d1_star,d2_star)
})

time2 = Sys.time()
# setwd("/panfs/jay/groups/4/ventz001/ghosh189/Opt_adap_design/New_Codes/results_debug")

setwd("/panfs/jay/groups/4/ventz001/ghosh189/Opt_adap_design/New_Codes/results2")

saveRDS(data_stage_12_full_list, paste0("NLP_dstar_new_mat_Sc",Job_mat[a,"Sc"],"batch",Job_mat[a,"Job"],"_c_",paste0(gsub("\\.", "point", as.character(constant[1]))),".rds"))
