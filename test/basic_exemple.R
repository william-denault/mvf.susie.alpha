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


indx_lst <- susiF.alpha::gen_wavelet_indx(lev_res = lev_res)
N = 100

#Number of covariates

P = 10

#Choosing which variable will have an effect
pos1 <- 1
pos2 <- 2




G = matrix(sample(c(0, 1,2), size=N*P, replace=T), nrow=N, ncol=P) #Genotype
beta1       <- 0
beta2       <- 0


noisy.data  <- list()
for ( i in 1:N)
{


  f1               <- effect1
  f2               <- effect2
  noisy.data [[i]] <-  t(beta1*G[i,pos1]*effect1 + beta2*G[i,pos2]*effect2 + matrix( rnorm((2^lev_res)*n_curve), ncol=n_curve))

}

noisy.data[[1]]
noisy.data[[2]]






dwt_data <- lapply(X = noisy.data, FUN = pack_dwt)

#line one in first matrix of dwt_data contain wt transform of condition 1 in ind 1
plot(dwt_data[[1]][1,], c(wd(noisy.data[[1]][1, ])$D,wd(noisy.data[[1]][1, ])$C[length(wd(noisy.data[[1]][1, ])$C)]))

plot(dwt_data[[1]][2,], c(wd(noisy.data[[1]][2, ])$D,
                          wd(noisy.data[[1]][2, ])$C[length(wd(noisy.data[[1]][2, ])$C)]))







DW_tens <- rearrange( dwt_data, lev_res = 7, n_curve=3)


dim(DW_tens)

X <- G


DW_tens[,1,] <- X[,1]+rnorm(N, sd=0.1)

tt <- mv_reg  (DW_tens[,1,],X[,1] )# correspond to regression of covariate 1 on a n dimensional wavelet coefficient
tt
bmv_reg(DW_tens[,1,],X[,1])


tt$Bhat/diag(tt$S)





### Fitting a mash per level of resolution

#idea generate two tensor, one coef one Shat
# then extract data to fit a mash

#need to parse by variable  (most inside loop) then by wave coef then by condition
v1 <- rep(1, N)


lol <- array(c( 10,2,0,0, 20,0,2,0,
         0,0,0,0), dim= c(2,2,3))





Y <- DW_tens

marg_assoc <- cal_Bhat_Shat_tensor  (Y, X, v1)

marg_assoc
hist(marg_assoc $tens_Bhat[,,1]/marg_assoc $tens_Shat[,,1], nclass = 100)
image(marg_assoc$tens_Bhat[,,1]/marg_assoc $tens_Shat[,,1])



fit_mash_level <- function(tens_marg, s, indx_lst)
{



   Bhat <- rbind_3Darray( tens_marg$tens_Bhat[,indx_lst[[s]],] )
   Shat <- rbind_3Darray( tens_marg$tens_Shat[,indx_lst[[s]],] )


   m <- basic_mash_fit(Bhat, Shat)

}
do.call( lapply(seq(dim(array)[3]), function(x)array[ , , x]))


Y <- DW_tens

Bhat  <- list()
Shat  <- list()

for ( l in 1:dim(Y)[2])
{
  out <-  do.call( cbind,lapply( 1:dim(X)[2], function(j) fit_lm(l= l,j=j, Y=Y,X=X, v1=v1  ) ) )
  Bhat[[l]]  <- out[1,]
  Shat[[l]] <- out[2,]
}

Bhat <- (do.call(cbind, Bhat))
Shat <- (do.call(cbind, Shat))
out  <- list( Bhat = Bhat,
              Shat = Shat)

lappl

Tens_coef # TxPx condition tensor


t(DW_tens[ ,1,])%*%X[,1]/(as.numeric(crossprod(X[,1])))


Y <- DW_tens[,1,]

dim(X)
dim(Y[ ,1,])

Bhat <- tt$Bhat
S <- tt$S
U <- diag( 1, 3)
S_inv <- solve (S)
multivariate_regression   (Bhat, S, U, S_inv)
