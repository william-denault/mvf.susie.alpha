#rm(list = ls())
library(ashr)
library(mashr)
library(abind)
library(wavethresh)
library(mvtnorm)
library(mixsqp)
library(mvf.susie.alpha)
set.seed(12)
n_curve=3
lev_res=5
effect1 <- mvf_susie_per_level(lev_res=lev_res,n_curve=3)$sim_func
effect2 <- mvf_susie_per_level(lev_res=lev_res,n_curve=3)$sim_func
verbose=TRUE
data.driven=FALSE
all= TRUE
indx_lst <- susiF.alpha::gen_wavelet_indx(lev_res = lev_res)
N = 100

#Number of covariates

P = 10

#Choosing which variable will have an effect
pos1 <- 1
pos2 <- 2




G = matrix(sample(c(0, 1,2), size=N*P, replace=T), nrow=N, ncol=P) #Genotype
beta1       <- 1
beta2       <- 1


noisy.data  <- list()
for ( i in 1:N)
{


  f1               <- effect1
  f2               <- effect2
  noisy.data [[i]] <-  t(beta1*G[i,pos1]*effect1 + beta2*G[i,pos2]*effect2 + matrix( rnorm(n= ((2^lev_res)*n_curve)), ncol=n_curve))

}

noisy.data[[1]]
noisy.data[[2]]
Y <- noisy.data

library(mvf.susie.alpha)



library(profvis)

profvis({
  out <- mvfsusie(Y=noisy.data,
                  X=X,
                  L=2, maxit = 10)
})
