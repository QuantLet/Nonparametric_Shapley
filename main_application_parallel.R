#BEGIN
bootstrap = function(j, B=100){

source("functions.R") # load all functions


library(np)
library(pracma)
library(cubature)
library(simstudy)
library(MASS)
#Get model fits and sort them in a list

sub_bw = subs; c=5
#d=2
sub_bw[[1]][1,1] = 4*c
sub_bw[[1]][1,2] = 38*c

sub_bw[[2]][1,1] = 4*c
sub_bw[[2]][2,1] = 40*c

model_list <<- model_list_fct(subs=subs, alt=FALSE, sub_bw) 



#generate B bootstrap samples
bs_shap <<- bs_shap=matrix(0, ncol=2, nrow=B)

#take bandwidths
b1 <<- model_list[[1]]$bw
b2 <<- model_list[[2]]$bw
b12 <<- model_list[[3]]$bw

dat_bs <<- cbind(Y,X)


ci_points = seq(120, 340, length.out=100)
ci_points2 = seq(2500, 5000, length.out=100)

for ( i in 1:B){
  print(i)
  model_list_bs = vector(mode="list")
  samp = dat[sample(nrow(dat_bs), B, replace=TRUE), ]
  model.np1 = npreg(reformulate("X1", "Y"), regtype = "ll", data = samp, bws=b1)
  model.np2 = npreg(reformulate("X2", "Y"), regtype = "ll", data = samp, bws=b2) 
  model.np12 = npreg(reformulate(c("X1", "X2"), "Y"), regtype = "ll", data = samp, bws=b12 ) 
  model_list_bs[[1]] = model.np1 
  model_list_bs[[2]] = model.np2 
  model_list_bs[[3]] = model.np12 
  
  bs_shap[i,1] = shapley_vec(j=1, c(ci_points[j], 2500), alt=FALSE, model_list = model_list_bs) 
  bs_shap[i,2] = shapley_vec(j=2, c(200, ci_points2[j]), alt=FALSE, model_list = model_list_bs) 
  
}


}
#closed


source("functions.R") # load all functions

#load data
f <- file.choose() # path noch festlegen!!

# Comma as separator and dot as decimal point by default
data=read.csv(file=f,                 # File name or full path of the file
              header = TRUE,        # Whether to read the header or not
              sep = ",",            # Separator of the values
              dec = ".")                  

data_imp = data[, c("price.baseMSRP","BLUETOOTH","HEATED.SEATS", "NAVIGATION", 
                    "features.Measurements.Curb.weight.lbs", "features.Engine.Horsepower.hp", "year", 
                    "features.Measurements.Width.in.", "features.Measurements.Length.in.", 
                    "features.Measurements.Height.in.", "Total.Seating", "features.Tires.and.Wheels.Tire.Aspect.Ratio")]

data_imp[,1] = data_imp[,1]/1000
names(data_imp) = c("price", "BLUETOOTH","heated","navi","weight", "hp", "year",
                    "width", "length", "height", "seats", "ratio")

#stratify
dat = data_imp[,c("price", "hp", "year", "weight")]

summary(dat)

dat_sort = dat[order(dat[,"year"]),]
dat_first = dat_sort[1:12230, ]
#dat_second = dat_sort[12231:24185, ]
#dat_third = dat_sort[24186:38435, ]

#All possible subsets
X <<- within(dat_first, rm("year", "price"))
d <<- ncol(X)
names = names(X)
names(X) = c("X1", "X2") 
Y <<- dat_first$price
subs <<- subsets(X)



library(parallel)
qs = mclapply(1:100, bootstrap, mc.cores=4)

points=matrix(unlist(qs), byrow = FALSE, ncol=100)


#extract
alpha = 0.1
qs = apply(bs_shap, 2, quantile, probs=c(alpha/2, 1 - alpha/2))

