#rm(list = ls())
library(ashr)
library(mashr)
library(tensorr)
library(wavethresh)
library(mvtnorm)
set.seed(1)
n_curve=3
lev_res=7
effect1 <- mvf_susie_per_level(lev_res=7,n_curve=3)$sim_func
effect2 <- mvf_susie_per_level(lev_res=7,n_curve=3)$sim_func

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
  noisy.data [[i]] <-  t(beta1*G[i,pos1]*effect1 + beta2*G[i,pos2]*effect2) #+ matrix( rnorm((2^lev_res)*n_curve), ncol=n_curve)

}

noisy.data[[1]]
noisy.data[[2]]




pack_dwt <- function( Y)
{
  W <- DWT2(Y)
  return(cbind( W$D,W$C))
}

dwt_data <- lapply(X = noisy.data, FUN = pack_dwt)

#line one in first matrix of dwt_data contain wt transform of condition 1 in ind 1
plot(dwt_data[[1]][1,], c(wd(noisy.data[[1]][1, ])$D,wd(noisy.data[[1]][1, ])$C[length(wd(noisy.data[[1]][1, ])$C)]))

plot(dwt_data[[1]][2,], c(wd(noisy.data[[1]][2, ])$D,
                          wd(noisy.data[[1]][2, ])$C[length(wd(noisy.data[[1]][2, ])$C)]))







DW_tens <- rearrange( dwt_data, lev_res = 7, n_curve=3)


dim(DW_tens)

X <- G


DW_tens[,1,] <- X[,1]+rnorm(N, sd=0.01)

mv_reg  (DW_tens[,1,],X[,1] )# correspond to regression of covariate 1 on a n dimensional wavelet coefficient

bmv_reg(DW_tens[,1,],X[,1])

t(DW_tens[ ,1,])%*%X[,1]/(as.numeric(crossprod(X[,1])))


Y <- DW_tens[,1,]

dim(X)
dim(Y[ ,1,])

bhat <- tt$bhat
S <- tt$S
U <- diag( 1, 3)
S_inv <- solve (S)
multivariate_regression   (bhat, S, U, S_inv)
