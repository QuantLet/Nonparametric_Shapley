library(parallel)
  
bs=function(m){
  ####Main code###
  #Load functions:
  library(np)
  library(pracma)
  library(cubature)
  library(simstudy)
  library(MASS)
  
  
  setwd("/Users/ratmir/MISE/2d_sim")
  
  
  source("functions.R") # load all functions
  
  cov=0
  sigma_sim=matrix(c(4, cov,
                     cov, 4), nrow=2, ncol=2)
  
  
  source("integral_population_2d.R")
  source("integral_estimation_2d.R") 
  
  source("shapley_int.R")
  #source("SE_vec_int.R")
  
  
  #separate functions for DGP. (additive case)
  g1 = function(X){ return( -sin(2*X[,1]) ) } #E(g1)=0. 
  g2 = function(X){ return( cos(2*X[,2])  ) } #E(g2)=0 
  #int = function(X){
  #  x1 = X[,1]
  #  x2 = X[,2]
  #  return( 2*cos(x1)*sin(2*x2)  ) 
  #}
  
  
  
  
  m_full_why = function(X){
    x1=as.numeric(X[1])
    x2=as.numeric(X[2])
    x3=as.numeric(X[3])
    return(-sin(2*x1) + cos(2*x2) )
  }
  
  
  true_model_list = list()
  true_model_list[[1]] = m_x1
  true_model_list[[2]] = m_x2
  true_model_list[[3]] = m_full_why 
  
  
  l = -2; u = 2; N=1000; M=1
  l_int = l; u_int = u
  
  d = 2
  
  M = 100
  #qs = matrix(0, ncol=2, nrow=M)
  
  
X<<-data.frame(mvrnorm(n=N, mu=c(0,0), Sigma=sigma_sim[1:2,1:2]))
#DGP

Y <<- g1(X) + g2(X)  + rnorm(nrow(X))
dat = data.frame(Y,X)
samp = vector(mode="list")

#generate B bootstrap samples
B=500
bs_shap=numeric()

#initial estimation
model.np.init1 = npreg(reformulate("X1", "Y"), regtype = "ll", data = dat, nmulti=0, bwmethod="cv.aic")
model.np.init2 = npreg(reformulate("X2", "Y"), regtype = "ll", data = dat, nmulti=0, bwmethod="cv.aic")
model.np.init12 = npreg(reformulate(c("X1", "X2"), "Y"), regtype = "ll", data = dat, nmulti=0, bwmethod="cv.aic")

#shapley estimator
model_list_init = vector(mode="list")
model_list_init[[1]] = model.np.init1 
model_list_init[[2]] = model.np.init2
model_list_init[[3]] = model.np.init12 

shap_estimated = shapley_vec(j=1, c(0,0), model_list = model_list_init) 

#take bandwidths
b1 = model.np.init1$bw
b2 = model.np.init2$bw
b12 = model.np.init12$bw

for ( i in 1:B){
  
model_list_bs = vector(mode="list")
  
samp = dat[sample(nrow(dat), 1000, replace=TRUE), ]


model.np1 = npreg(reformulate("X1", "Y"), regtype = "ll", data = samp, bws=b1)
model.np2 = npreg(reformulate("X2", "Y"), regtype = "ll", data = samp, bws=b2) 
model.np12 = npreg(reformulate(c("X1", "X2"), "Y"), regtype = "ll", data = samp, bws=b12 ) 

model_list_bs[[1]] = model.np1 
model_list_bs[[2]] = model.np2 #da
model_list_bs[[3]] = model.np12 #da

bs_shap[i] = shapley_vec(j=1, c(0,0), model_list = model_list_bs) #component-based 


}

#quantile approach
# take alpha/2 and 1 - alpha/2 quantile
#alpha = 0.05
#qs = quantile(bs_shap, prob=c(alpha/2, 1 - alpha/2))
sd_bs = sd(bs_shap)

# SE approach
CI_lower = shap_estimated - qnorm(1 - 0.05 / 2)*sd_bs
CI_upper = shap_estimated + qnorm(1 - 0.05 / 2)*sd_bs
qs = c(CI_lower, CI_upper)



return(qs)
}



library(parallel)
results_list = mclapply(1:100, bs, mc.cores=4)
results = matrix(unlist(results_list), byrow = FALSE, ncol=100)
qs = t(results)

#check if true shapley value is in this quantile range (better after MC runs)
shap_true = shapley_popul(j=1, c(0,0), model_list = model_list_init) #-> definiere model_list_init global
sum(shap_true >= qs[,1] & shap_true <= qs[,2])/M # larger equal than lower quantile, smaller equal than upper quantile
















#Y <<- g1(X) + g2(X)  + int(X) + rt(n=nrow(X), df=10)
#All possible subsets
#subs <<- subsets(X)
#Get model fits and sort them in a list, before:
#model_list <<- model_list_fct(subs=subs) #4sek für N=1000
dat = data.frame(Y,X)

















#Make a 3D plot to double check that shapley curve is correctly calculated.
# go in shapley_popul_vec(j, x_eval), where x_eval is 


x1_grid = seq(-2, 2, length.out=30) 
x2_grid = seq(-2, 2, length.out=30)


grid=t(expand.grid(x1_grid, x2_grid))
#evaluate on pop

shap_eval_est1 = matrix(0, nrow=ncol(grid), ncol=1)
shap_eval_est2 = matrix(0, nrow=ncol(grid), ncol=1)

for (i in 1:ncol(grid)){
  print(i)
  grid_col=as.numeric(grid[,i])
 shap_eval_est1[i] = shapley_int(j=1, grid_col) 
 shap_eval_est2[i] = shapley_int(j=2, grid_col) 
 
}

#shap_eval1=shapley_vec(j=1, grid) 
#shap_eval2=shapley_vec(j=2, grid)

#evaluate on est
shap_eval_est1=shapley_vec(j=1, grid) 
shap_eval_est2=shapley_vec(j=2, grid)

#shap_eval1=shapley_vec(j=1, grid) 
#shap_eval2=shapley_vec(j=2, grid)


shap_eval1 = matrix(0, nrow=ncol(grid), ncol=1)
shap_eval2 = matrix(0, nrow=ncol(grid), ncol=1)

for (i in 1:ncol(grid)){
  print(i)
  grid_col=as.numeric(grid[,i])
  shap_eval1[i] = shapley_popul(j=1, grid_col) 
  shap_eval2[i] = shapley_popul(j=2, grid_col) 
  
}

#shap_eval_est1=shapley_popul_vec(j=1, grid) 
#shap_eval_est2=shapley_popul_vec(j=2, grid)


shap_SE1 = SE_vec(j=1, grid)
surface1_SE=t(pracma::Reshape(shap_SE1, length(x1_grid), length(x2_grid))) #strong guess:columns are varied over x1 

shap_SE2 = SE_vec(j=2, grid)
surface2_SE=t(pracma::Reshape(shap_SE2, length(x1_grid), length(x2_grid))) #strong guess:columns are varied over x1 


surface1=t(pracma::Reshape(shap_eval1, length(x1_grid), length(x2_grid))) #strong guess:columns are varied over x1 
surface2=t(pracma::Reshape(shap_eval2, length(x1_grid), length(x2_grid))) #strong guess:columns are varied over x1 
surface1_est=t(pracma::Reshape(shap_eval_est1, length(x1_grid), length(x2_grid))) #strong guess:columns are varied over x1 
surface2_est=t(pracma::Reshape(shap_eval_est2, length(x1_grid), length(x2_grid))) #strong guess:columns are varied over x1 

surface1_SE = (surface1 - surface1_est)^2
surface2_SE = (surface2 - surface2_est)^2


library(plotly)
par(mar = c(0, 0, 0, 0))
fig1 = plot_ly(x=sort(x1_grid), y=sort(x2_grid), z=surface1, type="surface") %>% hide_colorbar()
fig1 <- fig1 %>% layout(scene=list(xaxis=list(title='x1'),yaxis=list(title='x2'),zaxis=list(title='shap1')))

fig1_est = plot_ly(x=sort(x1_grid), y=sort(x2_grid), z=surface1_est, type="surface") %>% hide_colorbar()
fig1_est <- fig1_est %>% layout(scene=list(xaxis=list(title='x1'),yaxis=list(title='x2'),zaxis=list(title='shap1')))

fig2 = plot_ly(x=sort(x1_grid), y=sort(x2_grid), z=surface2, type="surface") %>% hide_colorbar()
fig2 <- fig2 %>% layout(scene=list(xaxis=list(title='x1'),yaxis=list(title='x2'),zaxis=list(title='shap2')))

fig2_est = plot_ly(x=sort(x1_grid), y=sort(x2_grid), z=surface2_est, type="surface") %>% hide_colorbar()
fig2_est <- fig2 %>% layout(scene=list(xaxis=list(title='x1'),yaxis=list(title='x2'),zaxis=list(title='shap2')))

fig1_SE = plot_ly(x=sort(x1_grid), y=sort(x2_grid), z=surface1_SE, type="surface") %>% hide_colorbar()
fig1_SE <- fig1_SE %>% layout(scene=list(xaxis=list(title='x1'),yaxis=list(title='x2'),zaxis=list(title='SE shap1')))

fig2_SE = plot_ly(x=sort(x1_grid), y=sort(x2_grid), z=surface2_SE, type="surface") %>% hide_colorbar()
fig2_SE <- fig2_SE %>% layout(scene=list(xaxis=list(title='x1'),yaxis=list(title='x2'),zaxis=list(title='SE shap2')))



plot_ly(showscale = FALSE) %>%
  add_surface( x=~sort(x1_grid), y=~sort(x2_grid), z = ~surface1 , opacity = 1, colorscale = list(c(0,1),c("rgb(0,3,160)","rgb(0,3,160)"))) %>%
  add_surface(x=~sort(x1_grid), y=~sort(x2_grid),z = ~surface1_est, opacity = 0.3, colorscale = list(c(0,1),c("rgb(255,107,184)","rgb(128,0,64)"))) %>%
  layout(scene=list(xaxis=list(title='x1'), yaxis=list(title='x2'), zaxis=list(title=''),
                    camera = list(eye = list(x = 1.25, y = 1.25, z = 1.75))    ))

plot_ly(showscale = FALSE) %>%
  add_surface( x=~sort(x1_grid), y=~sort(x2_grid), z = ~surface2 , opacity = 1, colorscale = list(c(0,1),c("rgb(0,3,160)","rgb(0,3,160)"))) %>%
  add_surface(x=~sort(x1_grid), y=~sort(x2_grid),z = ~surface2_est, opacity = 0.3, colorscale = list(c(0,1),c("rgb(255,107,184)","rgb(128,0,64)"))) %>%
  layout(scene=list(xaxis=list(title='x1'), yaxis=list(title='x2'), zaxis=list(title=''),
                    camera = list(eye = list(x = 1.25, y = 1.25, z = 1.75))      ))



#SE

plot_ly(showscale = FALSE) %>%
  add_surface( x=~sort(x1_grid), y=~sort(x2_grid), z = ~surface1_SE , opacity = 1, colorscale = list(c(0,1),c("rgb(0,3,160)","rgb(0,3,160)"))) %>%
  layout(scene=list(xaxis=list(title='x1'), yaxis=list(title='x2'), zaxis=list(title='', range=c(0,0.42)),
                    camera = list(eye = list(x = 1.25, y = 1.25, z = 1.75))    )) # range=c(0,0.4) für zaxis range

plot_ly(showscale = FALSE) %>%
  add_surface( x=~sort(x1_grid), y=~sort(x2_grid), z = ~surface2_SE , opacity = 1, colorscale = list(c(0,1),c("rgb(0,3,160)","rgb(0,3,160)"))) %>%
  layout(scene=list(xaxis=list(title='x1'), yaxis=list(title='x2'), zaxis=list(title='', range=c(0,0.42)),
                    camera = list(eye = list(x = 1.25, y = 1.25, z = 1.75))    ))



fig1_SE
fig2_SE






