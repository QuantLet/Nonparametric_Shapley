library(parallel)


ISE_fct = function(m, obs){
  ####Main code###
  #Load functions:
  library(np)
  library(pracma)
  library(cubature)
  library(simstudy)
  library(MASS)
  
  source("functions.R") # load all functions
  
  #define cov before reading the functions
  cova<<-0
  sigma_sim<<-matrix(c(4, cova, cova,
                   cova, 4, cova,
                   cova, cova, 4), nrow=3, ncol=3)
  #Component-based:
  #CASE1: Independence
  source("integral_population.R")
  source("integral_estimation.R")
  source("shapley_int.R")
  source("SE_vec_int.R")
  
  
  #separate functions for DGP.
  g1 <<- function(X){ return( -sin(2*X[,1]) ) } #E(g1)=0. 
  g2 <<- function(X){ return( cos(3*X[,2])   ) } #E(g2)=0 
  g3 <<- function(X){ return( 0.5*X[,3] ) } #E(g3)=0. 
  int = function(X){
    x1 = X[,1]
    x2 = X[,2]
    return( 2*cos(x1)*sin(2*x2) ) 
  }
  
  
#  m_full_why <<- function(X){
#    x1=as.numeric(X[1])
#    x2=as.numeric(X[2])
#    x3=as.numeric(X[3])
#    return(-sin(2*x1) + cos(x2) + x3)
#  }
  
 # true_model_list = list()
#  true_model_list[[1]] = m_x1
 # true_model_list[[2]] = m_x2
  #true_model_list[[3]] = m_x3 
  #true_model_list[[4]] = m_x1_x2
  #true_model_list[[5]] = m_x1_x3
  #true_model_list[[6]] = m_x2_x3 
  #true_model_list[[7]] = m_full_why 
  
  l <<- -2; u <<- 2; N<<- obs
  l_int <<- l; u_int <<- u
  d <<- 3
  
  
  X<<-data.frame(mvrnorm(n=N, mu=c(0,0,0), Sigma=sigma_sim))
  
  #DGP
  Y <<- g1(X) + g2(X) + g3(X) + int(X) + rt(n=nrow(X), df=5)
  
  #All possible subsets
  subs <<- subsets(X)
  
  #Get model fits and sort them in a list
  model_list <<- model_list_fct(subs=subs) # 75sek, 14 sek with tol = 0.1 and ftol = 0.1 and 1 multistart
  
  while (sum(model_list[[7]]$bw[1:2]>10)>0
         | sum(model_list[[6]]$bw[1]>10)>0 
         | sum(model_list[[5]]$bw[1]>10)>0
         | sum(model_list[[4]]$bw[1:2]>10)>0
         | sum(model_list[[2]]$bw[1]>10)>0
         | sum(model_list[[1]]$bw[1]>10)>0)
  {
    X<<-data.frame(mvrnorm(n=N, mu=c(0,0,0), Sigma=sigma_sim))
    Y <<- g1(X) + g2(X) + g3(X) + int(X) + rt(n=nrow(X), df=5)
    model_list <<- model_list_fct(subs=subs)
  }
  
  
  
  # Component-based
  ISE_res1=hcubature(f=SE_vec, rep(l_int, d), rep(u_int, d), tol=3e-1, j=1)
  ISE_res2=hcubature(f=SE_vec, rep(l_int, d), rep(u_int, d), tol=3e-1, j=2)
  ISE_res3=hcubature(f=SE_vec, rep(l_int, d), rep(u_int, d), tol=3e-1, j=3)
  
  ISE1 = ISE_res1$integral 
  ISE2 = ISE_res2$integral 
  ISE3 = ISE_res3$integral 
  

  
  
  
  # Integral-based
  ISE_res1_int=hcubature(f=SE_vec_int, rep(l_int, d), rep(u_int, d), tol=3e-1, j=1)# 45sek für 0.3, gleich 0.5, intern zusätzl 0.3 gibt 31
  ISE_res2_int=hcubature(f=SE_vec_int, rep(l_int, d), rep(u_int, d), tol=3e-1, j=2)# 47sek für 0.3, gleich 0.5, intern zusätzl 0.3 gibt 34
  ISE_res3_int=hcubature(f=SE_vec_int, rep(l_int, d), rep(u_int, d), tol=3e-1, j=3)# 100sek für 0.3, 15 sek 0.5, intern zusätzl 0.3 gibt 11
  
  ISE1_int = ISE_res1_int$integral 
  ISE2_int = ISE_res2_int$integral  
  ISE3_int = ISE_res3_int$integral
  
  
  #print(c(ISE1_int, ISE2_int, ISE3_int))
  
  
  return(c(ISE1, ISE1_int, ISE2, ISE2_int, ISE3, ISE3_int))
}

res1=mclapply(1:200, ISE_fct, mc.cores=40, obs=300)
res2=mclapply(1:200, ISE_fct, mc.cores=40, obs=500)
res3=mclapply(1:200, ISE_fct, mc.cores=40, obs=1000)
res4=mclapply(1:200, ISE_fct, mc.cores=40, obs=2000)


results1=matrix(unlist(res1), byrow = FALSE, ncol=200)
results2=matrix(unlist(res2), byrow = FALSE, ncol=200)
results3=matrix(unlist(res3), byrow = FALSE, ncol=200)
results4=matrix(unlist(res4), byrow = FALSE, ncol=200)

a=rowMeans(results1)
b=rowMeans(results2)
c=rowMeans(results3)
d=rowMeans(results4)

tab_final=rbind(a,b,c,d)
