#rm(list = ls())
library(ashr)
library(mashr)
library(abind)
library(wavethresh)
library(mvtnorm)
library(mixsqp)
library(susiF.alpha)
library(mvf.susie.alpha)
set.seed(12)

min_levres <- 4

#gene expression effect of SNP 1 and SNP2
f01  <-  1
f02  <- -1


#effect of SNP 1 and SNP2 on mark 1
lev_res1 <- 5
temp_func <-  simu_IBSS_per_level(lev_res1)
f11 <-  temp_func$sim_func
temp_func <-  simu_IBSS_per_level(lev_res1)
f12 <-  temp_func$sim_func
#effect of SNP 1 and SNP2 on mark 2
lev_res2 <- 7
temp_func <-  simu_IBSS_per_level(lev_res2)
f21 <-  temp_func$sim_func
temp_func <-  simu_IBSS_per_level(lev_res2)
f22 <-  temp_func$sim_func

#effect of SNP 1 and SNP2 on 2 CpGs

f31  <-  c(1,-1)
f32  <-  c(-1,1)





indx_lst1 <- susiF.alpha::gen_wavelet_indx(lev_res = lev_res1)
indx_lst2 <- susiF.alpha::gen_wavelet_indx(lev_res = lev_res2)

N = 100

#Number of covariates

P = 10

#Choosing which variable will have an effect
pos1 <- 1




G = matrix(sample(c(0, 1,2), size=N*P, replace=T), nrow=N, ncol=P) #Genotype
beta1       <- 1



noisy.data  <- list()
for ( i in 1:N)
{

  Y0               <- f01*G[i,pos1]+f02*G[i,pos2] +rnorm(1) #gene expression
  Y1               <- t( beta1*G[i,pos1]*f11+ beta1*G[i,pos2]*f12+rnorm(2^lev_res1))#mark 1
  Y2               <- t( beta1*G[i,pos1]*f21+ beta1*G[i,pos2]*f22 +rnorm(2^lev_res2))#mark 2
  Y3               <-  f31*G[i, pos1]+f32*G[i, pos2] + c(rnorm(2))
  noisy.data [[i]] <-  list(Y0,Y1,Y2,Y3)

}

noisy.data[[1]]
noisy.data[[2]]
for ( k in 1:length(noisy.data[[1]]))
{
  list_dfs [[k]]     <- do.call(rbind, lapply(noisy.data, `[[`, k))
}

type_mark <-  is.functional(list_dfs)
type_mark
list_dfs  <- list()
list_wdfs <- list()
list_indx_lst  <-  list()
for ( k in which(type_mark=="functional"))
{
  list_dfs [[k]]     <- do.call(rbind, lapply(noisy.data, `[[`, k))
  temp               <- DWT2(list_dfs[[k]])
  list_wdfs[[k]]     <- cbind( temp$D,temp$C)
  list_indx_lst[[k]] <- gen_wavelet_indx( log2(ncol(  list_wdfs[[k]]) ))

}
Y_u <-
Y_f <- list_wdfs
v1 <- nrow(Y[[1]])
X <-G


G_prior <- init_prior_multfsusie(Y,
                                 X,
                                 v1,
                                 list_indx_lst
                                 )
multf


dwt_data <- lapply(X = noisy.data, FUN = pack_dwt)

#line one in first matrix of dwt_data contain wt transform of condition 1 in ind 1
plot(dwt_data[[1]][1,], c(wd(noisy.data[[1]][1, ])$D,wd(noisy.data[[1]][1, ])$C[length(wd(noisy.data[[1]][1, ])$C)]))

plot(dwt_data[[1]][2,], c(wd(noisy.data[[1]][2, ])$D,
                          wd(noisy.data[[1]][2, ])$C[length(wd(noisy.data[[1]][2, ])$C)]))







DW_tens <- rearrange( dwt_data, lev_res = lev_res, n_curve=3)


dim(DW_tens)

X <- G

