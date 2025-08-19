rm(list=ls())
library(parallel)
library(MASS)
setwd("/panfs/jay/groups/4/ventz001/ghosh189/Opt_adap_design/New_Codes")
source("imp_functions.R")
setwd("/panfs/jay/groups/4/ventz001/ghosh189/Opt_adap_design/New_Codes/results")


N_sims=1000
N1 = 100
N2 = 100

N  =  N1+N2

Scenarios=c(1:17)
futil_thresh = 0.2 ###################### similar to other futlity stopping criteria


mclapply(Scenarios,mc.cores=4,function(Sc){
  set.seed(Sc) 
  
  time1=Sys.time()
  
  prop     = scenarios_prop_phi_nlp(Sc)$prop
  b.phi      = scenarios_prop_phi_nlp(Sc)$b.phi
  
  
  bx345 = c(0.5,0.5,0.5) 
  b.phi = c(b.phi,bx345)
  phi        = b.phi[3:6]
  
  prob_T_Rct =2/3
  
  rho1     =1
  lambda      = 5
  
  indices     = c(1:4)
  n.subsets   = length(indices)
  
  #phi         = rep(0,n.subsets)
  subsets     = lapply(1:(2^n.subsets - 1), function(i) indices[which(bitwAnd(i, 2^(0:(n.subsets - 1))) != 0)])
  
  #phi       = c(phi.1=-0.1, phi.2=-0.3, phi.3=1, phi.4=1.3)
  weighted.phi.H0.true         = sapply(subsets,function(indices.vec)(prop[indices.vec]) %*% phi[indices.vec] /sum(prop[indices.vec])) 
  size.subsets <- lapply(subsets,function(indices.vec) sum(prop[indices.vec])^rho1)
  
  
  v.Y       = 1
  
  
  #' ============ @data-stage-1
  
  n_stage1_trials = 1000
  
  
  data_stage1_list=replicate(n = n_stage1_trials , simplify=F, expr = {
    
    indices= c(1,2,3,4)
    R.i       = data_gen2(indices,N,prop,prob_T_Rct)
    Y.i.H0    =  R.i%*%b.phi + rnorm(N1)
    
    list(Y.i.H0=Y.i.H0, R.i=R.i)
  })
  
  
  # pvalues_list = lapply(1:n_stage1_trials, function(i){
  #   
  #   D1=data_stage1_list[[i]]
  #   # posterior from stage 1
  #   
  #   R.i       = D1$R.i
  #   Y.i.H0    = D1$Y.i.H0
  #   
  #   R.t.R     = t(R.i) %*% R.i
  #   R.t.Y.H0  = t(R.i) %*% Y.i.H0
  #   
  #   var     = solve( R.t.R)      
  #   mu.hat   = var %*% (R.t.Y.H0)
  #   
  #   TE.H0                       = c(mu.hat[3],mu.hat[4],mu.hat[5],mu.hat[6])
  #   TE.var                        =var[3:6,3:6]
  #   # Number of observations (n) and dimension of TE.H0 (p)
  #   n = nrow(R.i)
  #   p = length(TE.H0)
  #   
  #   # Hotelling's T-squared test statistic
  #   T2 = t(TE.H0) %*% solve(TE.var) %*% TE.H0
  #   T2 = as.numeric(T2)
  #   
  #   Hotel_tsq_stat  = T2
  #   df = p
  #   p.value         =pchisq(Hotel_tsq_stat, df, lower.tail = FALSE)
  #   # Convert T-squared statistic to an F-statistic
  #   F.stat = (n - p) / (p * (n - 1)) * T2
  # 
  #   # Degrees of freedom for the F-distribution
  #   df1 = p
  #   df2 = n - p
  # 
  #   # p-value from the F-distribution
  #   p.value = pf(F.stat, df1, df2, lower.tail = FALSE)
  # 
  #   list(p.value)
  # })
  
  
  pvalues_list = sapply(1:n_stage1_trials, function(i){
    
    D1=data_stage1_list[[i]]
    # posterior from stage 1
    
    R.i       = D1$R.i[1:N1,]
    Y.i.H0    = D1$Y.i.H0[1:N1]
    
    trt_indic = apply(R.i[,3:6]==1,1,any)
    
    #model <- lm(Y.i.H0 ~ R.i[, c(1, 2, 7, 8, 9)] + trt_indic)
    model <- lm(Y.i.H0 ~ trt_indic)
    coef_summary = summary(model)$coefficients
    trt_coef = coef_summary["trt_indicTRUE", "Estimate"]
    trt_se = coef_summary["trt_indicTRUE", "Std. Error"]
    
    # Calculating tail probability that P(coefficient of TE > 0)
    prob = 1 - pnorm(0, mean = trt_coef, sd = trt_se)
    
    if(prob<futil_thresh){
      d1=16
      p_value_ztest =1
      p_value_lm=1
    }else{
      d1=15
    R.i       = D1$R.i
    Y.i.H0    = D1$Y.i.H0
    
    ctrl_indic = apply(R.i[,3:6]==0,1,all)
    z_test    = t.test(Y.i.H0[!ctrl_indic],Y.i.H0[ctrl_indic],alternative = "greater",var.equal = TRUE)
    p_value_ztest     = z_test$p.value
    model= summary(lm(Y.i.H0~as.numeric(1-ctrl_indic)+R.i[,c(1:2,7:9)]-1))
    zstat =  model$coefficients[1,1]/model$coefficients[1,2]
    p_value_lm = 1-pnorm(zstat)
    
    }
    c(d1,p_value_ztest,p_value_lm)
  })
  
  not_early_stop = which(pvalues_list[1,]!=16)
 
  pvals = rowMeans(pvalues_list[2:3,not_early_stop]<0.05)
  my_list   = list(power_ztest=pvals[1],power_lm=pvals[2])
  
  stage_1_decisions= pvalues_list[1,]
  
  ### ztest#########################################################
  final_d2_stars = as.numeric(pvalues_list[2,not_early_stop]<0.05)
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
  saveRDS(my_list,paste0("RCT_ztest",Sc,"_new.rds")   )
  
  
  
  ### lm#########################################################
  final_d2_stars = as.numeric(pvalues_list[3,not_early_stop]<0.05)
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
  saveRDS(my_list,paste0("RCT_lm",Sc,"_new.rds")   )
})



