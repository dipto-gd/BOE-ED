rm(list=ls())


.libPaths("/home/faird/ghosh189/R/x86_64-pc-linux-gnu-library/4.3/")
library("StanHeaders")
library("rstan")

library(MASS)
library(dplyr)
library(mvtnorm,lib.loc="/home/faird/ghosh189/R/x86_64-pc-linux-gnu-library/4.0")
setwd("/projects/standard/ventz001/ghosh189/Opt_adap_design/New_Codes")
source("imp_functions.R")

no_scenarios = 17


N_trials = 1000
no_jobs = 50

n_stage1_trials = round(N_trials/no_jobs)
n_stage1_trials=1
Job_mat = expand.grid(Job=1:no_jobs,Sc=1:(no_scenarios))

# job runs from no_jobs *no_scenarios
a = 3
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


  
  #'===========@internal_data
  indices.vec= c(1,2,3,4)
  
  p3  = p4 = p5 = 0.5
  X.i.2.labels    = sample(as.character(indices.vec), N1,replace = TRUE,prob =  prop[indices.vec])
  X.i.2.labels    = as.numeric(X.i.2.labels)
  
  values          = c(0,1)
  X.i.2           =  expand.grid(rep(list(values), 2))
  X.i_list      =  lapply(X.i.2.labels,function(i)X.i.2[i,])
  X.i_mat       =  do.call(rbind,X.i_list)
  
  
  # Convert Z.i rows to G.i rows
  G.i.2       = t(apply(X.i_mat, 1, convert_to_G))
  T.i.2       = rbinom(N1,1,prob_T_Rct)
  L.i.2       = (G.i.2  * T.i.2)
  X.i.3       = cbind(rbinom(N1,1,p3),rbinom(N1,1,p4),rbinom(N1,1,p5))
  
  R.i.2       = as.matrix(cbind(X.i_mat,  L.i.2,X.i.3))
  
  R.i_RCT       = R.i.2
  
  Y.i_RCT    =  R.i_RCT%*%b.phi + rnorm(N1)
  
  
  RCT_data = data.frame( Y= Y.i_RCT, 
                         X1= X.i_mat[,1],
                         X2= X.i_mat[,2],
                           X3= X.i.3[,1],
                           X4=X.i.3[,2],
                           X5= X.i.3[,3] )
  #prop of control is given by
  # mean(apply(R.i, 1, function(row) all(row[3:6] == 0)))
  
  #'===========@external_data
  indices= c(1,2,3,4)
  
  p3  = p4 = p5 = 0.5
  X.i.2.labels    = sample(as.character(indices.vec), N_ed,replace = TRUE,prob =  prop[indices.vec])
  X.i.2.labels    = as.numeric(X.i.2.labels)
  
  values          = c(0,1)
  X.i.2           =  expand.grid(rep(list(values), 2))
  X.i_list      =  lapply(X.i.2.labels,function(i)X.i.2[i,])
  X.i_mat       =  do.call(rbind,X.i_list)
  
  
  # Convert Z.i rows to G.i rows
  G.i.2       = t(apply(X.i_mat, 1, convert_to_G))
  T.i.2       = rbinom(N_ed,1,prob_T_Ed)
  L.i.2       = (G.i.2  * T.i.2)
  X.i.3       = cbind(rbinom(N_ed,1,p3),rbinom(N_ed,1,p4),rbinom(N_ed,1,p5))
  
  R.i.2       = as.matrix(cbind(X.i_mat,  L.i.2,X.i.3))
  R.i_ED       = data_gen2(indices,N_ed,prop_ED,prob_T_Ed)
  
  Y.i_ED    = intercept+ R.i_ED%*%b.phi_ED + rnorm(N_ed)
  
  
  ED_data = data.frame( Y= Y.i_ED, 
                         X1= X.i_mat[,1],
                         X2= X.i_mat[,2],
                         X3= X.i.3[,1],
                         X4=X.i.3[,2],
                         X5= X.i.3[,3] )

  save(RCT_data, file = "RCT_data.RData")
  save(ED_data, file = "ED_data.RData")
  