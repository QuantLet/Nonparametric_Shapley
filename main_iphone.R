####Main code###
#Load functions:
library(np)
library(pracma)
library(cubature)
library(simstudy)
library(MASS)


setwd("/Users/ratmir/MISE")


source("functions.R") # load all functions

#define cov before reading the functions
cov=0
sigma_sim=matrix(c(4, cov, cov,
                   cov, 4, cov,
                   cov, cov, 4), nrow=3, ncol=3)
#define model parameters before reading the functions
#beta_01 <<- 1; beta_02 <<- 1; beta_11 <<- 0.5; beta_12 <<- 0.5
#gamm <<- 5; con <<- 0

#Component-based:
#CASE1: Independence
source("integral_population.R")
source("integral_estimation.R")
source("shapley_int.R")
source("SE_vec_int.R")


#separate functions for DGP. (additive case)
g1 = function(X){ return( pmax(0, 1 - abs(X[,1])/0.2 ) ) } #E(g1)=0. 
g2 = function(X){ return( pmax(0, 1 - abs(X[,2])/0.6 )  ) } #E(g2)=0 
g3 = function(X){ return( pmax(0, 1 - abs(X[,3])/0.9 ) ) } #E(g3)=0. 
#int = function(X){
#  x1 = X[,1]
#  x2 = X[,2]
#  return( 2*cos(x1)*sin(2*x2)  ) 
#}



m_full_why = function(X){
  x1=as.numeric(X[1])
  x2=as.numeric(X[2])
  x3=as.numeric(X[3])
  return(pmax(0, 1 - abs(x1)/0.2 ) + pmax(0, 1 - abs(x2)/0.6 ) + pmax(0, 1 - abs(x3)/0.9 ) )
}


true_model_list = list()
true_model_list[[1]] = m_x1
true_model_list[[2]] = m_x2
true_model_list[[3]] = m_x3 
true_model_list[[4]] = m_x1_x2
true_model_list[[5]] = m_x1_x3
true_model_list[[6]] = m_x2_x3 
true_model_list[[7]] = m_full_why 



l = -2; u = 2; N=1000; M=1
l_int = l; u_int = u
d = 3




ISE_fct = function(m){
  
  X<<-data.frame(mvrnorm(n=N, mu=c(0,0,0), Sigma=sigma_sim))
  
  #DGP
  Y <<- g1(X) + g2(X) + g3(X) + rnorm(nrow(X), mean=0, sd=1)
  #Y <<- g1(X) + g2(X) + g3(X) + int(X) + rnorm(nrow(X), mean=0, sd=1)
  
  
  
  #All possible subsets
  subs <<- subsets(X)
  
  #Get model fits and sort them in a list
  model_list <<- model_list_fct(subs=subs) 
  
  while (sum(model_list[[7]]$bw[1:2]>10)>0
         | sum(model_list[[6]]$bw[1]>10)>0 
         | sum(model_list[[5]]$bw[1]>10)>0
         | sum(model_list[[4]]$bw[1:2]>10)>0
         | sum(model_list[[2]]$bw[1]>10)>0
         | sum(model_list[[1]]$bw[1]>10)>0)
  {
    X<<-data.frame(mvrnorm(n=N, mu=c(0,0,0), Sigma=sigma_sim))
    Y <<- g1(X) + g2(X) + g3(X) + rnorm(nrow(X))
    model_list <<- model_list_fct(subs=subs)
  }
  
  
  
  
  
  
  print(model_list[[7]])
  # Component-based
  ISE_res1=hcubature(f=SE_vec, rep(l_int, d), rep(u_int, d), tol=3e-1, j=1)
  ISE_res2=hcubature(f=SE_vec, rep(l_int, d), rep(u_int, d), tol=3e-1, j=2)
  ISE_res3=hcubature(f=SE_vec, rep(l_int, d), rep(u_int, d), tol=3e-1, j=3)
  
  ISE1 = ISE_res1$integral 
  ISE2 = ISE_res2$integral 
  ISE3 = ISE_res3$integral 
  
  #print(c(ISE1, ISE2, ISE3))
  
  
  #system.time(hcubature(f=SE_vec_int, rep(l_int, d), rep(u_int, d), tol=3e-1, j=1))[3]
  #system.time(hcubature(f=SE_vec_int, rep(l_int, d), rep(u_int, d), tol=1e-1, j=1))[3]
  
  #system.time(hcubature(f=SE_vec_int, rep(l_int, d), rep(u_int, d), tol=3e-1, j=2))[3]
  #system.time(hcubature(f=SE_vec_int, rep(l_int, d), rep(u_int, d), tol=1e-1, j=2))[3]
  
  #system.time(hcubature(f=SE_vec_int, rep(l_int, d), rep(u_int, d), tol=3e-1, j=3))[3]
  #system.time(hcubature(f=SE_vec_int, rep(l_int, d), rep(u_int, d), tol=1e-1, j=3))[3]
  
  
  
  
  # Integral-based
  ISE_res1_int=hcubature(f=SE_vec_int, rep(l_int, d), rep(u_int, d), tol=3e-1, j=1)# 45sek für 0.3, gleich 0.5, intern zusätzl 0.3 gibt 31
  ISE_res2_int=hcubature(f=SE_vec_int, rep(l_int, d), rep(u_int, d), tol=3e-1, j=2)# 47sek für 0.3, gleich 0.5, intern zusätzl 0.3 gibt 34
  ISE_res3_int=hcubature(f=SE_vec_int, rep(l_int, d), rep(u_int, d), tol=3e-1, j=3)# 100sek für 0.3, 15 sek 0.5, intern zusätzl 0.3 gibt 11
  
  ISE1_int = ISE_res1_int$integral 
  ISE2_int = ISE_res2_int$integral  
  ISE3_int = ISE_res3_int$integral
  
  
  #print(c(ISE1_int, ISE2_int, ISE3_int))
  
  return(c(ISE1, ISE2, ISE1_int, ISE2_int))
  #return(1)
  
}



x1_grid = seq(-1, 1, length.out=200) 
x2_grid = seq(-1, 1, length.out=200)


grid=t(expand.grid(x1_grid, x2_grid))
grid = rbind(grid, rep(0,ncol(grid)))


shap_eval1 = matrix(0, nrow=ncol(grid), ncol=1)
shap_eval2 = matrix(0, nrow=ncol(grid), ncol=1)
#shap_eval3 = matrix(0, nrow=ncol(grid), ncol=1)

for (i in 1:ncol(grid)){
  print(i)
  grid_col=as.numeric(grid[,i])
  shap_eval1[i] = shapley_popul(j=1, grid_col) 
  shap_eval2[i] = shapley_popul(j=2, grid_col) 
  #shap_eval3[i] = shapley_popul(j=3, grid_col) 
  
}




surface1=t(pracma::Reshape(shap_eval1, length(x1_grid), length(x2_grid))) #strong guess:columns are varied over x1 
surface2=t(pracma::Reshape(shap_eval2, length(x1_grid), length(x2_grid))) #strong guess:columns are varied over x1 
#surface3=t(pracma::Reshape(shap_eval3, length(x1_grid), length(x2_grid))) #strong guess:columns are varied over x1 


library(plotly)
par(mar = c(0, 0, 0, 0))
fig1 = plot_ly(x=sort(x1_grid), y=sort(x2_grid), z=surface1, type="surface") %>% hide_colorbar()
fig1 <- fig1 %>% layout(scene=list(xaxis=list(title='x1'),yaxis=list(title='x2'),zaxis=list(title='', range=c(-0.01, 1.1) ), 
                                   camera = list(eye = list(x = 1.25, y = 1.25, z = 1.75))) )
fig1

fig2 = plot_ly(x=sort(x1_grid), y=sort(x2_grid), z=surface2, type="surface") %>% hide_colorbar()
fig2 <- fig2 %>% layout(scene=list(xaxis=list(title='x1'),yaxis=list(title='x2'),zaxis=list(title='', range=c(-0.01, 1.1) ),
                                   camera = list(eye = list(x = 1.25, y = 1.25, z = 1.75))))
fig2

fig3 = plot_ly(x=sort(x1_grid), y=sort(x2_grid), z=surface3, type="surface") %>% hide_colorbar()
fig3 <- fig3 %>% layout(scene=list(xaxis=list(title='x1'),yaxis=list(title='x3'),zaxis=list(title='shap3')))






