setwd("/Users/ratmir/MISE/application")
source("functions.R") # load all functions

cov=0
sigma_sim=matrix(c(1, cov,
                   cov, 1), nrow=2, ncol=2)

source("integral_estimation_2d.R") 
source("shapley_int.R")

library(np)
library(pracma)
library(cubature)
library(simstudy)
library(MASS)
#load data
f <- file.choose()

# Comma as separator and dot as decimal point by default
data=read.csv(file=f,                 # File name or full path of the file
              header = TRUE,        # Whether to read the header or not
              sep = ",",            # Separator of the values
              dec = ".")                  

data_imp = data[, c("price.baseMSRP","BLUETOOTH","HEATED.SEATS", "NAVIGATION", 
"features.Measurements.Curb.weight.lbs", "features.Engine.Horsepower.hp", "year", 
"features.Measurements.Width.in.", "features.Measurements.Length.in.", 
"features.Measurements.Height.in.", "Total.Seating", "features.Tires.and.Wheels.Tire.Aspect.Ratio", 
"features.Measurements.Wheel.base.in.", "features.Tires.and.Wheels.Wheel.Diameter")]
data_imp[,1] = data_imp[,1]/1000

#standardize the cont. variables
#data_imp[,3] = (data_imp[,3] - mean(data_imp[,3])) / sd(data_imp[,3])
#data_imp[,4] = (data_imp[,4] - mean(data_imp[,4])) / sd(data_imp[,4])
#data_imp[,5] = (data_imp[,5] - mean(data_imp[,5])) / sd(data_imp[,5])
#data_imp[,6] = (data_imp[,6] - mean(data_imp[,6])) / sd(data_imp[,6])
#data_imp[,7] = (data_imp[,7] - mean(data_imp[,7])) / sd(data_imp[,7])
#data_imp[,8] = (data_imp[,8] - mean(data_imp[,8])) / sd(data_imp[,8])
#data_imp[,9] = (data_imp[,9] - mean(data_imp[,9])) / sd(data_imp[,9])


#bianry
data_imp$BLUETOOTH = ifelse(data_imp$BLUETOOTH == 'Yes', 1, 0)
data_imp$BLUETOOTH = as.factor(data_imp$BLUETOOTH)
data_imp$HEATED.SEATS = ifelse(data_imp$HEATED.SEATS == 'Yes', 1, 0)
data_imp$HEATED.SEATS = as.factor(data_imp$HEATED.SEATS)
data_imp$NAVIGATION = ifelse(data_imp$NAVIGATION == 'Yes', 1, 0)
data_imp$NAVIGATION = as.factor(data_imp$NAVIGATION)

#categorical
#data_imp$CAR.TYPE = as.factor(data_imp$CAR.TYPE)
#levels(data_imp$CAR.TYPE) = c(1:6)

names(data_imp) = c("price", "BLUETOOTH","heated","navi","weight", "hp", "year",
                    "width", "length", "height", "seats", "ratio", "wheels", "diameter")


#stratify
dat = data_imp[,c("price", "hp", "year", "weight", "length")]

summary(dat)

dat_sort = dat[order(dat[,"year"]),]
#dat_sort = dat_sort[!(dat_sort[,2]<120 | dat_sort[,2] >400),]
#dat_sort = dat_sort[!(dat_sort[,4]<1990 | dat_sort[,4] >6000),]
#dat_sort = dat_sort[!(dat_sort[,5]<170 | dat_sort[,5] >260),]





dat_first = dat_sort[1:12230, ]
dat_second = dat_sort[12231:24185, ]
dat_third = dat_sort[24186:38435, ]

#split by year
a=sum(dat_sort[order(dat_sort[,"year"]),"year"]==2001)
b=sum(dat_sort[order(dat_sort[,"year"]),"year"]==2002)
c=sum(dat_sort[order(dat_sort[,"year"]),"year"]==2003)
d=sum(dat_sort[order(dat_sort[,"year"]),"year"]==2004)
e=sum(dat_sort[order(dat_sort[,"year"]),"year"]==2005)
f=sum(dat_sort[order(dat_sort[,"year"]),"year"]==2006)
g=sum(dat_sort[order(dat_sort[,"year"]),"year"]==2007)
#12230
aa=sum(dat_sort[order(dat_sort[,"year"]),"year"]==2008)
bb=sum(dat_sort[order(dat_sort[,"year"]),"year"]==2009)
cc=sum(dat_sort[order(dat_sort[,"year"]),"year"]==2010)
dd=sum(dat_sort[order(dat_sort[,"year"]),"year"]==2011)
ee=sum(dat_sort[order(dat_sort[,"year"]),"year"]==2012)
ff=sum(dat_sort[order(dat_sort[,"year"]),"year"]==2013)
#11955
aaa=sum(dat_sort[order(dat_sort[,"year"]),"year"]==2014)
bbb=sum(dat_sort[order(dat_sort[,"year"]),"year"]==2015)
ccc=sum(dat_sort[order(dat_sort[,"year"]),"year"]==2016)
ddd=sum(dat_sort[order(dat_sort[,"year"]),"year"]==2017)
eee=sum(dat_sort[order(dat_sort[,"year"]),"year"]==2018)
fff=sum(dat_sort[order(dat_sort[,"year"]),"year"]==2019)
ggg=sum(dat_sort[order(dat_sort[,"year"]),"year"]==2020)
#14250

#All possible subsets
X = within(dat_first, rm("year", "price"))

#X = within(dat_sort, rm("year", "price"))

d = ncol(X)
names = names(X)
names(X) = c("X1", "X2", "X3") 
Y = dat_first$price
#Y = dat_sort$price

#summary statistics
par(mar = c(1.3, 2.5, 3, 0.4))
boxplot(X[!(X[,1]>400 | X[,1]<120),1], cex.main=3,
        cex.axis=2, font.main = 2, main="Horsepower")
boxplot(X[!(X[,2]<2000 | X[,2]>6000),2], cex.main=3,
        cex.axis=2, font.main = 2, main="Weight")
boxplot(X[!(X[,3]<170 | X[,3]>260),3], cex.main=3,
        cex.axis=2, font.main = 2, main="Length")








subs <<- subsets(X)

#Get model fits and sort them in a list

sub_bw = subs; c=1
#d=2
#sub_bw[[1]][1,1] = 4*c
#sub_bw[[1]][1,2] = 38*c

#sub_bw[[2]][1,1] = 4*c
#sub_bw[[2]][2,1] = 40*c

#d=3 für diameter
#sub_bw[[1]][,1] = 5*c
#sub_bw[[1]][,2] = 55*c
#sub_bw[[1]][,3] = 0.2*c

#sub_bw[[2]][1,1] = 5*c
#sub_bw[[2]][2,1] = 55*c

#sub_bw[[2]][1,2] = 5*c
#sub_bw[[2]][2,2] = 0.2*c

#sub_bw[[2]][1,3] = 55*c
#sub_bw[[2]][2,3] = 0.2*c

#sub_bw[[3]][1,1] = 5*c
#sub_bw[[3]][2,1] = 55*c
#sub_bw[[3]][3,1] = 0.2*c

c=3  #current for length
sub_bw[[1]][,1] = 1*c
sub_bw[[1]][,2] = 40*c
sub_bw[[1]][,3] = 2*c # ok

sub_bw[[2]][1,1] = 4*c
sub_bw[[2]][2,1] = 60*c

sub_bw[[2]][1,2] = 24*c
sub_bw[[2]][2,2] = 6*c

sub_bw[[2]][1,3] = 80*c
sub_bw[[2]][2,3] = 1*c # ?

sub_bw[[3]][1,1] = 12*c
sub_bw[[3]][2,1] = 80*c
sub_bw[[3]][3,1] = 4*c










model_list <<- model_list_fct(subs=subs, alt=FALSE, sub_bw) 





#grid_a = seq(summary(X[,1])[1], summary(X[,1])[6], length.out=100)
grid_a = seq(120, 400, length.out=100)
grid_b = seq(2000, 6000, length.out=100)
grid_c = seq(170, 250, length.out=100)

#grid_c = seq(summary(X[,3])[1], 240, length.out=100)
#grid_d = seq(summary(60, 90, length.out=100))



eval_a = rbind(grid_a, rep(3500, 100),  rep(190, 100))
col_a = shapley_vec(j=1, x_eval=eval_a , alt=FALSE, model_list = model_list)

eval_b = rbind(rep(250, 100), grid_b, rep(190, 100))
col_b = shapley_vec(j=2, x_eval=eval_b , alt=FALSE, model_list = model_list)

eval_c = rbind(rep(250, 100), rep(3500, 100), grid_c) #200/250 PS
col_c = shapley_vec(j=3, x_eval=eval_c, alt=FALSE, model_list = model_list)

#col_a_int = numeric()
#col_b_int = numeric()

#for(i in 1:ncol(eval_a) ){
#  print(i)
#col_a_int[i] = shapley_int(j=1, x_eval = eval_a[,i])
#col_b_int[i] = shapley_int(j=2, x_eval = eval_b[,i])
#
#}




#eval_c = rbind(rep(mean_a, 100), rep(mean_b, 100), grid_c, rep(mean_d, 100))
#col_c = shapley_vec(j=3, x_eval=eval_c , alt=TRUE)


#Delete outliers
#outl_a = col_a<=100 & col_a>=-100
#clean_grid_a = grid_a[which(outl_a)]
#clean_col_a = col_a[outl_a]

#outl_b = col_b<=70 & col_b>=-70
#clean_grid_b = grid_b[which(outl_b)]
#clean_col_b = col_b[outl_b]

#outl_c = col_c<=70 & col_c>=-70
#clean_grid_c = grid_c[which(outl_c)]
#clean_col_c = col_c[outl_c]

par(mar = c(2.5, 2.5, 2.2, 0.4))
plot(y=col_a, x=grid_a, type="l", xlab="", ylab="", ylim=c(-13,70), main="2001 - 2007", cex.main=2,
     cex.axis=1.5, font.main = 2)
x1_ci = seq(120, 400, length.out=48)
points(x=x1_ci, y=points[1,] ,cex=0.3, col="red")
points(x=x1_ci, y=points[2,] ,cex=0.3, col="red")
lines(x=x1_ci, y=points[1,] ,cex=0.3, col="red")
lines(x=x1_ci, y=points[2,] ,cex=0.3, col="red")
abline(a=0,b=0, lty=2)


points(x=x1_ci, y=points[7,] ,cex=0.3, col="blue")
points(x=x1_ci, y=points[8,] ,cex=0.3, col="blue")
lines(x=x1_ci, y=points[7,] ,cex=0.3, col="blue")
lines(x=x1_ci, y=points[8,] ,cex=0.3, col="blue")

plot(y=col_b, x=grid_b, type="l", xlab="X2", ylab="", ylim=c(-11,22), main="2008 - 2013", cex.main=2,
     cex.axis=1.5, font.main = 2)
x2_ci = seq(2000, 6000, length.out=48)
points(x=x2_ci, y=points[3,] ,cex=0.3, col="red")
points(x=x2_ci, y=points[4,] ,cex=0.3, col="red")
lines(x=x2_ci, y=points[3,] ,cex=0.3, col="red")
lines(x=x2_ci, y=points[4,] ,cex=0.3, col="red")
abline(a=0,b=0, lty=2)

points(x=x2_ci, y=points[9,] ,cex=0.3, col="blue")
points(x=x2_ci, y=points[10,] ,cex=0.3, col="blue")
lines(x=x2_ci, y=points[9,] ,cex=0.3, col="blue")
lines(x=x2_ci, y=points[10,] ,cex=0.3, col="blue")

plot(y=col_c, x=grid_c, type="l", xlab="X3", ylab="", ylim=c(-25,15), main="2008 - 2013", cex.main=2,
     cex.axis=1.5, font.main = 2)
x3_ci = seq(170, 260, length.out=48)
points(x=x3_ci, y=points[5,] ,cex=0.3, col="red")
points(x=x3_ci, y=points[6,] ,cex=0.3, col="red")
lines(x=x3_ci, y=points[5,] ,cex=0.3, col="red")
lines(x=x3_ci, y=points[6,] ,cex=0.3, col="red")

abline(a=0,b=0, lty=2)


points(x=x3_ci, y=points[11,] ,cex=0.3, col="blue")
points(x=x3_ci, y=points[12,] ,cex=0.3, col="blue")
lines(x=x3_ci, y=points[11,] ,cex=0.3, col="blue")
lines(x=x3_ci, y=points[12,] ,cex=0.3, col="blue")
#3D plots

x1_grid = seq(150, 330, length.out=100) 
x2_grid = seq(2000, 6000, length.out=100)


grid = t(expand.grid(x1_grid, x2_grid))
grid = rbind(grid , rep(190,ncol(grid)) )

shap_eval_est1=shapley_vec(j=1, grid, alt=FALSE, model_list = model_list) 
shap_eval_est2=shapley_vec(j=2, grid, alt=FALSE, model_list = model_list)

#On data points
#shap_eval_est1=shapley_vec(j=1, t(X), alt=TRUE, model_list = model_list) 
#shap_eval_est2=shapley_vec(j=2, t(X), alt=TRUE, model_list = model_list)

#plot(y=shap_eval_est1, x=X[,1])
#plot(y=shap_eval_est2, x=X[,2])



#exclude outliers

outl_one = shap_eval_est1<=100 & shap_eval_est1>=-100
outl_two = shap_eval_est2<=100 & shap_eval_est2>=-100


shap_eval_est1[!outl_one] = median(shap_eval_est1)
shap_eval_est2[!outl_two] = median(shap_eval_est2)





#On X3:

x2_grid = seq(2000, 6000, length.out=100) 
x3_grid = seq(170, 250, length.out=100)

#grid = t(expand.grid(x1_grid, x3_grid))
#grid = rbind(rep(250,ncol(grid)),grid)

grid = t(expand.grid(x2_grid, x3_grid))
grid = rbind(rep(250,ncol(grid)), grid[1,], grid[2,])


shap_eval_est3=shapley_vec(j=3, grid, alt=FALSE, model_list = model_list) 


surface3_est=t(pracma::Reshape(shap_eval_est3, length(x1_grid), length(x3_grid))) 





surface1_est=t(pracma::Reshape(shap_eval_est1, length(x1_grid), length(x2_grid))) #strong guess:columns are varied over x1 
surface2_est=t(pracma::Reshape(shap_eval_est2, length(x1_grid), length(x2_grid))) #strong guess:columns are varied over x1 



library(plotly)
par(mar = c(0, 0, 0, 0))
fig1_est = plot_ly(x=sort(x1_grid), y=sort(x2_grid), z=surface1_est, type="surface") %>% hide_colorbar()
fig1_est <- fig1_est %>% layout(scene=list(xaxis=list(title='x1'),yaxis=list(title='x2'),zaxis=list(title='shap1')))

fig2_est = plot_ly(x=sort(x1_grid), y=sort(x2_grid), z=surface2_est, type="surface") %>% hide_colorbar()
fig2_est <- fig2_est %>% layout(scene=list(xaxis=list(title='x1'),yaxis=list(title='x2'),zaxis=list(title='shap2')))

fig3_est = plot_ly(x=sort(x2_grid), y=sort(x3_grid), z=surface3_est, type="surface") %>% hide_colorbar()
fig3_est <- fig1_est %>% layout(scene=list(xaxis=list(title='x1'),yaxis=list(title='x3'),zaxis=list(title='shap3')))


mrg <- list(l = 0, r = 0,
            b = 0, t = 0,
            pad = 0)

plot_ly(showscale = FALSE) %>%
  add_surface(x=~sort(x1_grid), y=~sort(x2_grid),z = ~surface1_est, opacity = 0.3, colorscale = list(c(0,1),c("rgb(255,107,184)","rgb(128,0,64)"))) %>%
  layout(scene=list(xaxis=list(title='Horsepower'), yaxis=list(title='Weight'), zaxis=list(title=''),
                    camera = list(eye = list(x = 1.25, y = 1.25, z = 1.75)) ), margin=mrg)
#different colors
plot_ly(showscale = FALSE) %>%
  add_surface(x=~sort(x1_grid), y=~sort(x2_grid),z = ~surface1_est, opacity = 0.8) %>%
  layout(scene=list(xaxis=list(title='Horsepower'), yaxis=list(title='Weight'), zaxis=list(title='', range = c(-20,50)),
                    camera = list(eye = list(x = 1.25, y = 1.25, z = 1.75)) ), margin=mrg)


plot_ly(showscale = FALSE) %>%
  add_surface(x=~sort(x1_grid), y=~sort(x2_grid),z = ~surface2_est, opacity = 0.3, colorscale = list(c(0,1),c("rgb(255,107,184)","rgb(128,0,64)"))) %>%
  layout(scene=list(xaxis=list(title='Horsepower'), yaxis=list(title='Weight'), zaxis=list(title=''),
                    camera = list(eye = list(x = 1.25, y = 1.25, z = 1.75)) ), margin=mrg)

#different colors
plot_ly(showscale = FALSE) %>%
  add_surface(x=~sort(x1_grid), y=~sort(x2_grid),z = ~surface2_est, opacity = 0.8) %>%
  layout(scene=list(xaxis=list(title='Horsepower'), yaxis=list(title='Weight'), zaxis=list(title='', range = c(-15,40)),
                    camera = list(eye = list(x = 1.25, y = 1.25, z = 1.75)) ), margin=mrg)


plot_ly(showscale = FALSE) %>%
  add_surface(x=~sort(x1_grid), y=~sort(x3_grid),z = ~surface3_est, opacity = 0.3, colorscale = list(c(0,1),c("rgb(255,107,184)","rgb(128,0,64)"))) %>%
  layout(scene=list(xaxis=list(title='Horsepower'), yaxis=list(title='Length'), zaxis=list(title=''),
                    camera = list(eye = list(x = 1.25, y = 1.25, z = 1.75)) ), margin=mrg)

#different colors
plot_ly(showscale = FALSE) %>%
  add_surface(x=~sort(x2_grid), y=~sort(x3_grid),z = ~surface3_est, opacity = 0.8) %>%
  layout(scene=list(xaxis=list(title='Horsepower'), yaxis=list(title='Length'), zaxis=list(title='', range = c(-20,8)),
                    camera = list(eye = list(x = 1.25, y = 1.25, z = 1.75)) ), margin=mrg)



#Points
plot_ly(type = "scatter3d", x = X[1:4000,1], y = X[1:4000,2], z = X[1:4000,3], mode = "markers")




#training set and test set

set.seed(101) # Set Seed so that same sample can be reproduced in future also
# Now Selecting 75% of data as sample from total 'n' rows of the data  
train <- sample(nrow(data_imp), 0.7*nrow(data_imp), replace = FALSE)

data_train = data_imp[train, ]
data_test  = data_imp[-train, ]

#y_train = data_train[,1]
#y_test = data_test[,1]


X_train = data_train[,-1]
X_test = data_test[,-1]

Y_train = data_train[,1]
Y_test = data_test[,1]

#All possible subsets
#subs <<- subsets(X)

#Get model fits and sort them in a list
#model_list <<- model_list_fct(subs=subs) 

#Diagnostic plots on training set.
plot(y=data_train$price, x=data_train$hp) # the higher the hp the higher the price variance
plot(y=data_train$price, x=data_train$weight) 
plot(y=data_train$price, x=data_train$year) # price doesnt increase that much over years
plot(y=data_train$price, x=data_train$BLUETOOTH) # bluetooth more expensive
plot(y=data_train$price, x=data_train$width) 
plot(y=data_train$price, x=data_train$length) 
plot(y=data_train$price, x=data_train$height) 
plot(y=data_train$price, x=data_train$seats) 
plot(y=data_train$price, x=data_train$type) 


#1. Start with OLS.
fit_l = lm(price ~ ., data=data_train)
#Residual plots: QQ plot, residual vs fitted. Can check for normality, linearity, homoskedasticity
plot(fit_l)


#How good is the model? Evaluate on independent, unseen test data set:
ols_mse_test = mean((predict(fit_l, newdata=X_test) - Y_test)^2)


# Indications for non-linearity and heteroscedasitcity. Let us use RF, since it is 
# more flexible.
library(randomForest)
model1 = randomForest(price ~ ., data = data_train, mtry=4, importance=TRUE)
rf_mse_test = mean((predict(model1, newdata = X_test) - Y_test)^2)
rf_mse_test
ols_mse_test 
#Plot fitted vs true Y on test to see how good the prediction is.
plot((predict(model1, newdata = X_test[1:200,])), x=Y_test[1:200], ylim=c(0,100), xlim=c(0,100))
abline(0,1)

#As we see, RF performs better than OLS on independent test set. Deploy the final model.























