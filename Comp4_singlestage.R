rm(list=ls())
library(parallel)
.libPaths("/home/faird/ghosh189/R/x86_64-pc-linux-gnu-library/4.2/")
library(MASS)
library(dplyr)
setwd("/panfs/jay/groups/4/ventz001/ghosh189/Opt_adap_design/New_Codes")
source("imp_functions.R")
setwd("/panfs/jay/groups/4/ventz001/ghosh189/Opt_adap_design/New_Codes/results")
#setwd("C:/Users/User/Desktop/Academic/UMN/Fall 2023/Steffen-RA/Simulations/Local Codes/local_results")



  a  = as.integer(commandArgs()[3])  # set to a=1, or a=2, ...., a=200

  time1=Sys.time()
  rho1        = 1
  rho2        = 0
  lambda      =5
  
  
  indices     = c(1:4)
  n.subsets   = length(indices)
  subsets     = lapply(1:(2^n.subsets - 1), function(i) indices[which(bitwAnd(i, 2^(0:(n.subsets - 1))) != 0)])
  subsets.string =  sapply(subsets, function(x) paste(x, collapse = ","))
  
  futil_thresh = 0.4 ###################### similar to other futlity stopping criteria
  
  Sc=a
   set.seed(Sc) 
  time1=Sys.time()
  
  prop     = scenarios_prop_phi_nlp(Sc)$prop
  prop_ED  = scenarios_prop_phi_nlp(Sc)$prop_ED
  b.phi      = scenarios_prop_phi_nlp(Sc)$b.phi
  b.phi_ED   = scenarios_prop_phi_nlp(Sc)$b.phi_ED
  
  
  bx345 = c(0.5,0.5,0.5) 
  b.phi = c(b.phi,bx345)
  b.phi_ED = c(b.phi_ED, bx345)
  
  phi        = b.phi[3:6]
  phi_ED     = b.phi_ED[3:6]
  
  prob_T_Rct =2/3
  
  size.subsets           = lapply(subsets,function(indices.vec) sum(prop[indices.vec])^rho1)
  
  weighted.phi.H0.true         = sapply(subsets,function(indices.vec)(prop[indices.vec]) %*% phi[indices.vec] /sum(prop[indices.vec])) 
  
  N1        = 100
  N2        = 100
  N_ed      = 500
  v.Y       = 1
  
  prior.Ivar = diag(1/100, length(b.phi))
  prior.mu   = rep(0,length(b.phi))
  
  
  #' ============ @data-stage-1
  time1=Sys.time()
  
  # Generate n_stage1_trials 
  n_stage1_trials = 1000
  
  
  
  
  N  = N1+N2
  #c  = c_vec[j]
  #N1 = as.integer(N/2)
  #N2 = N-N1
  data_stage1_list=replicate(n = n_stage1_trials , simplify=F, expr = {
    
    
    #' ====================== @internal_data
    # Generate all populations with prop proportions
    
    R.i       = data_gen2(indices,N1+N2,prop,prob_T_Rct)
    
    Y.i    =  R.i%*%b.phi + rnorm(N1+N2)
    
    
    list(Y.i.H0=Y.i, R.i=R.i)
  })
  
 
#' ============ @loop_over_stage1_datasets

#' following gives optimal decision at the end of stage 1
data_stage_12_full_list = lapply(1:n_stage1_trials, function(i){
  
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
    d1_star=16
    d2_star=0
  }else{
  

 
  
  # STAGE 2, USE FULL DATA
  R.i       = D1$R.i
  Y.i.H0    = D1$Y.i.H0
  R.t.R     = t(R.i) %*% R.i
  R.t.Y.H0  = t(R.i) %*% Y.i.H0
  post.var     = solve(prior.Ivar + R.t.R)      
  post.mu.H0   = post.var %*% (prior.Ivar %*% prior.mu + R.t.Y.H0)

  TE.2.H0.ext  = post.mu.H0[3:6]
  TE.var.2.ext = post.var[3:6,3:6]

  weighted.phi         = sapply(subsets,function(indices.vec)(prop[indices.vec]) %*% TE.2.H0.ext[indices.vec] /sum(prop[indices.vec])) 
  
  #' true treatment effects
  util.subsets.H0.true              = util_subsets_func(phi,subsets,prop,rho1,rho2)
  #' expected treatment effects
  TE.H0                       = c(post.mu.H0[3],post.mu.H0[4],post.mu.H0[5],post.mu.H0[6])
  util.subsets.H0.int            = util_subsets_func(TE.H0,subsets,prop,rho1,rho2)
  
  
  d1_star              = which.max(util.subsets.H0.int)
  
  indices.vec=subsets[[d1_star]]
  phi.j.mu.H0  = prop[indices.vec] %*% TE.2.H0.ext[indices.vec] /sum(prop[indices.vec]) 
  phi.j.var = prop[indices.vec] %*% TE.var.2.ext[indices.vec,indices.vec]%*% (prop[indices.vec]) /sum(prop[indices.vec])^2
  
  c_j_star    = -0.41
  d2_star     = ifelse(phi.j.mu.H0 >-c_j_star,1,0) 
  

  }

  list(d2_star=d2_star,d1_star=d1_star)
  
})


weighted.phi.H0.true=c(0,weighted.phi.H0.true)

stage_1_decisions           = sapply(data_stage_12_full_list, function(V)V$d1_star)

table(stage_1_decisions)

final_d2_stars               = sapply(data_stage_12_full_list, function(V)V$d2_star)

table(final_d2_stars)

#' prop of early stop is 0 since design does not stop early
#' 
#" What proportion of trials are stopped early?
prop_earlystop           =  mean(stage_1_decisions==16)

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
subpops_with_atleast_one_positive = which(sapply(subsets,function(V){ifelse(any(phi[V]>0),1,0)})==1)
stage_2_part_benefit        = which( stage_1_decisions %in% subpops_with_atleast_one_positive) 
stage_2_accept_part_benefit = which(final_d2_stars[stage_2_part_benefit]==1)
if(length(stage_2_part_benefit)>0){
  power    = length(stage_2_accept_part_benefit)/length(stage_1_decisions)
}else{power=0}

#' Proportion of times a drug is accepted when everyone in population is harmed (Type 1 error)
subpops_with_all_negative = which(sapply(subsets,function(V){ifelse(all(phi[V]<=0),1,0)})==1)
stage_2_all_harm        = which( stage_1_decisions %in% subpops_with_all_negative)
stage_2_accept_all_harm = which(final_d2_stars[which( stage_1_decisions %in% subpops_with_all_negative)]==1)

if(length(stage_2_all_harm)>0){
  type_1_error    = length(stage_2_accept_all_harm)/length(stage_1_decisions)
}else{
  type_1_error=0
}

#' Proportion of times a drug is accepted when no one in population is harmed
subpops_with_no_negative = which(sapply(subsets,function(V){ifelse(all(phi[V]>=0),1,0)})==1)
stage_2_no_harm        = which( stage_1_decisions %in% subpops_with_no_negative)
stage_2_accept_no_harm = which(final_d2_stars[which( stage_1_decisions %in% subpops_with_no_negative)]==1)
if(length(stage_2_no_harm)>0){
  prop_accept_no_harm    = length(stage_2_accept_no_harm)/length(stage_1_decisions)
}else{
  prop_accept_no_harm=0
}

#' Proportion of times a drug is accepted when everyone in population benefits
subpops_with_all_pos = which(sapply(subsets,function(V){ifelse(all(phi[V]>0),1,0)})==1)
stage_2_all_benft        = which( stage_1_decisions %in% subpops_with_all_pos)
stage_2_accept_all_benft = which(final_d2_stars[which( stage_1_decisions %in% subpops_with_all_pos)]==1)
if(length(stage_2_all_benft)>0){
  prop_accept_all_benft    = length(stage_2_accept_all_benft)/length(stage_1_decisions)
}else{
  prop_accept_all_benft=0
}

#' Proportion of subgroups with TE>0 discovered
subpops_with_atleast_one_pos = which(sapply(subsets,function(V){ifelse(any(phi[V]>0),1,0)})==1)
final_accept_index = which(final_d2_stars==1)
right_selection_prop = sapply(final_accept_index, function(dummy_index){
  ifelse(stage_1_decisions[dummy_index] %in% subpops_with_atleast_one_pos,
         sum(prop[which(subpops_with_atleast_one_pos %in% subsets[[stage_1_decisions[dummy_index]]])])/sum(prop[which(phi>0)]) ,0)
  
})
right_selection_prop = sum(right_selection_prop)/1000

#' Proportion of subgroups with TE<=0  discovered

subpops_with_atleast_one_neg = which(sapply(subsets,function(V){ifelse(any(phi[V]<=0),1,0)})==1)
final_accept_index = which(final_d2_stars==1)
wrong_selection_prop = sapply(final_accept_index, function(dummy_index){
  ifelse(stage_1_decisions[dummy_index] %in% subpops_with_atleast_one_neg,
         sum(prop[which(subpops_with_atleast_one_neg %in% subsets[[stage_1_decisions[dummy_index]]])])/sum(prop[which(phi<=0)]) ,0)
  
})
wrong_selection_prop = sum(wrong_selection_prop)/1000


# assume lambda patients arrive in the center per month.Add 2 weeks to get test results outcomes from arrival.
len_enroll  =mean(sapply(stage_1_decisions,function(V){sum(rexp(N1,rate=lambda))}))+mean(sapply(stage_1_decisions,function(V){ sum(rexp(N2,rate=lambda*(sum(prop[unlist(subsets[V])]))))}),na.rm=TRUE) + 0.5

my_list   = list(type_1_error=type_1_error,
                 power=power,
                 prop_accept_no_harm=prop_accept_no_harm,
                 prop_accept_all_benft=prop_accept_all_benft,
                 right_selection_prop=right_selection_prop,
                 wrong_selection_prop= wrong_selection_prop,
                 len_enroll =len_enroll,
                 prop_earlystop=prop_earlystop)

my_list
saveRDS(my_list,paste0("Comp4_Sc_",Sc,"N2_",N,".rds")   )




