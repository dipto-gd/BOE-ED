rm(list=ls())

setwd("/panfs/jay/groups/4/ventz001/ghosh189/Opt_adap_design/New_Codes")
source("imp_functions.R")
setwd("/panfs/jay/groups/4/ventz001/ghosh189/Opt_adap_design/New_Codes/results")
#setwd("C:/Users/User/Desktop/Academic/UMN/Fall 2023/Steffen-RA/Simulations/Local Codes/local_results")
library(parallel)
.libPaths("/home/faird/ghosh189/R/x86_64-pc-linux-gnu-library/4.2/")
library(MASS)

N1        = 100
N2        = 100
#N_ed      = 500
v.Y       = 1

N  = N1+N2


lapply(1:17,function(a){
  Sc=a
file_path <- paste0("Comp3_Sc",Sc,"samplesize_",N,"u0.rds")

if (file.exists(file_path)) {
  file.remove(file_path)
} 
})

a           = as.integer(commandArgs()[3])  # set to a=1, or a=2, ...., a=200

time1=Sys.time()
rho1        = 1
rho2        = 0
lambda      = 5

scale       = 1

indices     = c(1:4)
n.subsets   = length(indices)
subsets     = lapply(1:(2^n.subsets - 1), function(i) indices[which(bitwAnd(i, 2^(0:(n.subsets - 1))) != 0)])
subsets.string =  sapply(subsets, function(x) paste(x, collapse = ","))

futil_thresh = 0.4 ###################### similar to other futlity stopping criteria
Sc=a
#set.seed(Sc) 
time1=Sys.time()

prop     = scenarios_prop_phi_nlp(Sc)$prop
prop_ED  = scenarios_prop_phi_nlp(Sc)$prop_ED
b.phi = scenarios_prop_phi_nlp(Sc)$b.phi
b.phi_ED   = scenarios_prop_phi_nlp(Sc)$b.phi_ED

bx345 = c(0.5,0.5,0.5) 
b.phi = c(b.phi,bx345)
b.phi_ED = c(b.phi_ED, bx345)

phi        = b.phi[3:6]
phi_ED     = b.phi_ED[3:6]

prob_T_Rct =2/3

size.subsets           = lapply(subsets,function(indices.vec) sum(prop[indices.vec])^rho1)

weighted.phi.H0.true         = sapply(subsets,function(indices.vec)(prop[indices.vec]) %*% phi[indices.vec] /sum(prop[indices.vec])) 



prior.Ivar = diag(1/100, length(b.phi))
prior.mu   = rep(0,length(b.phi))


#' ============ @data-stage-1
time1=Sys.time()
EPPTE_EPNTE = EPPTE_EPNTE_Function_OPY()
# Generate n_stage1_trials 
n_stage1_trials = 1000

  #' ============ @loop_over_stage1_datasets_under_null_to_get_cutoff
  #' 
  # see when sc 17-20 occur what sould change to get the u0s for them????????????????????????????????????????????????????????????????/
  if (all(phi)==0 & !file.exists(paste0("Comp3_Sc",Sc,"samplesize_",N,"u0.rds"))){
 
    data_stage1_list=replicate(n = n_stage1_trials , simplify=F, expr = {
      
      
      #' ====================== @internal_data
      # Generate all populations with prop proportions
      
      R.i       = data_gen2(indices,N1+N2,prop,prob_T_Rct)
      
      Y.i    =  R.i%*%b.phi + rnorm(N1+N2)
      
      
      list(Y.i.H0=Y.i, R.i=R.i)
    })
   non_futile_utilities = sapply(1:n_stage1_trials, function(i){
      
      D1=data_stage1_list[[i]]
      
      R.i       = D1$R.i
      Y.i.H0    = D1$Y.i.H0
      R.t.R     = t(R.i) %*% R.i
      R.t.Y.H0  = t(R.i) %*% Y.i.H0
      post.var     = solve(prior.Ivar + R.t.R)      
      post.mu.H0   = post.var %*% (prior.Ivar %*% prior.mu + R.t.Y.H0)
      post.mu.H0
      
      #' true treatment effects
      util.subsets.H0.true              = util_subsets_func(phi,subsets,prop,rho1,rho2)
      #' expected treatment effects
      TE.H0                       = c(post.mu.H0[3],post.mu.H0[4],post.mu.H0[5],post.mu.H0[6])
      util.subsets.H0.int            = util_subsets_func(TE.H0,subsets,prop,rho1,rho2)
      
      util1_15max = max(util.subsets.H0.int)
      util1_15max
    })
    
    u0  = quantile(non_futile_utilities,0.95,lower.tail=TRUE)
    saveRDS(u0,paste0("Comp3_Sc",0,"samplesize_",N,"u0.rds"))
    
  }
    
  data_stage1_list=replicate(n = n_stage1_trials , simplify=F, expr = {
    
    
    #' ====================== @internal_data
    # Generate all populations with prop proportions
    
    R.i       = data_gen2(indices,N1+N2,prop,prob_T_Rct)
    
    Y.i    =  R.i%*%b.phi + rnorm(N1+N2)
    
    
    list(Y.i.H0=Y.i, R.i=R.i)
  })
    u0 =readRDS(paste0("Comp3_Sc",0,"samplesize_",N,"u0.rds"))/scale
    data_stage_12_full_list = sapply(1:n_stage1_trials, function(i){
      
      D1=data_stage1_list[[i]]
      
      R.i       = D1$R.i[1:N1,]
      Y.i.H0    = D1$Y.i.H0[1:N1,]
      R.t.R     = t(R.i) %*% R.i
      R.t.Y.H0  = t(R.i) %*% Y.i.H0
      post.var     = solve(prior.Ivar + R.t.R)      
      post.mu.H0   = post.var %*% (prior.Ivar %*% prior.mu + R.t.Y.H0)
      
      TE.2.H0.ext  = post.mu.H0[3:6]
      TE.var.2.ext = post.var[3:6,3:6]
      
      prob = 1- pnorm(0.05, mean= TE.2.H0.ext,sd= sqrt(diag(TE.var.2.ext)))
      
      if(all(prob<futil_thresh)){
        d1_star=17
        d2_star=0
      }else{
        
      
      R.i       = D1$R.i
      Y.i.H0    = D1$Y.i.H0
      R.t.R     = t(R.i) %*% R.i
      R.t.Y.H0  = t(R.i) %*% Y.i.H0
      post.var     = solve(prior.Ivar + R.t.R)      
      post.mu.H0   = post.var %*% (prior.Ivar %*% prior.mu + R.t.Y.H0)
      post.mu.H0
      
      #' true treatment effects
      util.subsets.H0.true              = util_subsets_func(phi,subsets,prop,rho1,rho2)
      #' expected treatment effects
      TE.H0                       = c(post.mu.H0[3],post.mu.H0[4],post.mu.H0[5],post.mu.H0[6])
      util.subsets.H0.int            = util_subsets_func(TE.H0,subsets,prop,rho1,rho2)
      
      util.subsets.H0.int = c(util.subsets.H0.int,u0)
      
      d1_star              = which.max(util.subsets.H0.int)
      
      }
      
      d1_star
      
      
    })
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  stage_1_decisions           = sapply(data_stage_12_full_list, function(V)V)
  
  table(stage_1_decisions)
  
  final_d2_stars               = sapply(data_stage_12_full_list, function(V) ifelse(V %in% c(16,17),0,1))
  
  prop_earlystop =  mean(stage_1_decisions==17)
  
  #' #'How many times is the drug correctly rejected? 
  #' subpops_with_leq0_effect     = which(weighted.phi.H0.true<=0)
  #' stage_2_avg_harm         = which( stage_1_decisions %in% subpops_with_leq0_effect)
  #' stage_2_avg_harm_accept  = which(unlist(final_d2_stars)[stage_2_avg_harm]==1)
  #' if(length(stage_2_avg_harm)>0){
  #'   type_1_error= length(stage_2_avg_harm_accept)/length(stage_1_decisions)
  #' }else{
  #'   type_1_error=0  
  #'   print("No pop with average treatment harm selected")
  #' }
  #' 
  #' #' How many times is the drug correctly accepted? 
  #' subpops_with_greater0_effect = which(weighted.phi.H0.true>0)
  #' stage_2_avg_benefit         = which( stage_1_decisions %in% subpops_with_greater0_effect)
  #' stage_2_avg_benefit_accept  = which(unlist(final_d2_stars)[stage_2_avg_benefit]==1)
  #' if(length(stage_2_avg_benefit)>0){
  #'   power= length(stage_2_avg_benefit_accept)/length(stage_1_decisions)
  #' }else{power=0}
  #' In whom is the drug correctly accepted?
  # stage_2_correct_accept_decisions = stage_1_decisions[which(final_d2_stars[which( stage_1_decisions %in% subpops_with_greater0_effect)]==1)]
  # 
  # table(stage_2_correct_accept_decisions)
  
  #' Proportion of times a drug is accepted when some parts of the subpopulation benefit (Power)
  subpops_with_atleast_one_positive = which(sapply(1:length(subsets),function(V){ifelse(any(phi_function(V,Sc)<0),1,0)})==1)
  # stage_2_part_benefit        = which( stage_1_decisions %in% subpops_with_atleast_one_positive) 
  # stage_2_accept_part_benefit = which(final_d2_stars[stage_2_part_benefit]==1)
  # if(length(stage_2_part_benefit)>0){
  #   power    = length(stage_2_accept_part_benefit)/length(stage_1_decisions)
  # }else{power=0}
  # 
  
  power = mean((final_d2_stars==1) & (stage_1_decisions %in% subpops_with_atleast_one_positive ))
  #' Proportion of times a drug is accepted when everyone in population is harmed (Type 1 error)
  subpops_with_all_negative = which(sapply(1:length(subsets),function(V){ifelse(all(phi_function(V,Sc)>=0),1,0)})==1)
  # stage_2_all_harm        = which( stage_1_decisions %in% subpops_with_all_negative)
  # stage_2_accept_all_harm = which(final_d2_stars[which( stage_1_decisions %in% subpops_with_all_negative)]==1)
  # 
  # if(length(stage_2_all_harm)>0){
  #   type_1_error    = length(stage_2_accept_all_harm)/length(stage_1_decisions)
  # }else{
  #   type_1_error=0
  #  }
  
  stage_2_all_harm        =  (stage_1_decisions %in% subpops_with_all_negative) & (final_d2_stars==1)
  type_1_error= mean(stage_2_all_harm)
  
  #' Proportion of times a drug is accepted when no one in population is harmed
  subpops_with_no_negative = which(sapply(subsets,function(V){ifelse(all(phi_function(V,Sc)<=0),1,0)})==1)
  # stage_2_no_harm        = which( stage_1_decisions %in% subpops_with_no_negative)
  # stage_2_accept_no_harm = which(final_d2_stars[which( stage_1_decisions %in% subpops_with_no_negative)]==1)
  # if(length(stage_2_no_harm)>0){
  #   prop_accept_no_harm    = length(stage_2_accept_no_harm)/length(stage_1_decisions)
  # }else{
  #   prop_accept_no_harm=0
  # }
  prop_accept_no_harm =  mean((final_d2_stars==1) & (stage_1_decisions %in% subpops_with_no_negative ))
  #' Proportion of times a drug is accepted when everyone in population benefits
  subpops_with_all_pos = which(sapply(subsets,function(V){ifelse(all(phi_function(V,Sc)<0),1,0)})==1)
  # stage_2_all_benft        = which( stage_1_decisions %in% subpops_with_all_pos)
  # stage_2_accept_all_benft = which(final_d2_stars[which( stage_1_decisions %in% subpops_with_all_pos)]==1)
  # if(length(stage_2_all_benft)>0){
  #   prop_accept_all_benft    = length(stage_2_accept_all_benft)/length(stage_1_decisions)
  # }else{
  #   prop_accept_all_benft=0
  # }
  prop_accept_all_benft = mean((stage_1_decisions %in% subpops_with_all_pos) & final_d2_stars==1)
  
  
  Sc_mod_5= ifelse(Sc>5,Sc-5,Sc) 
  #' Proportion of subgroups with TE>0 discovered
  PTE = EPPTE_EPNTE[[Sc_mod_5]][1,]
  subpops_with_atleast_one_pos = which(sapply(1:(length(subsets)-1),function(V){ifelse(any(PTE[V]>0),1,0)})==1)
  final_accept_index = which(final_d2_stars==1)
  right_selection_prop = sapply(final_accept_index, function(dummy_index){
    ifelse(stage_1_decisions[dummy_index] %in% subpops_with_atleast_one_pos,1 ,0)
    
  })
  right_selection_prop = sum(right_selection_prop)/1000
  
  #' Proportion of subgroups with TE<=0  discovered
  NTE = EPPTE_EPNTE[[Sc_mod_5]][2,]
  subpops_with_atleast_one_neg = which(sapply(1:(length(subsets)-1),function(V){ifelse(any(NTE[V]<=0),1,0)})==1)
  final_accept_index = which(final_d2_stars==1)
  wrong_selection_prop = sapply(final_accept_index, function(dummy_index){
    ifelse(stage_1_decisions[dummy_index] %in% subpops_with_atleast_one_neg, 1,0)
    
  })
  wrong_selection_prop = sum(wrong_selection_prop)/1000
  
  
  # assume lambda patients arrive in the center per month.Add 2 weeks to get test results outcomes from arrival.
  len_enroll  =mean(sapply(stage_1_decisions,function(V){sum(rexp(N1,rate=lambda))}))+mean(sapply(stage_1_decisions,function(V){ sum(rexp(N2,rate=lambda*(subsets_prev[V])))}),na.rm=TRUE) + 0.5
  my_list   = list(type_1_error=type_1_error,
                   power=power,
                   prop_accept_no_harm=prop_accept_no_harm,
                   prop_accept_all_benft=prop_accept_all_benft,
                   # right_selection_prop=right_selection_prop,
                   # wrong_selection_prop= wrong_selection_prop,
                   # len_enroll =len_enroll,
                   prop_earlystop=prop_earlystop)
  
  NLP_noED_list[[Sc]] = my_list 
  print(Sc)
  
  
  my_list
  saveRDS(my_list,paste0("Comp3_Sc_",Sc,"N2_",N,".rds")   )
 

