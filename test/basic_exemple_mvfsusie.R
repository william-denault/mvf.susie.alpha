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






dwt_data <- lapply(X = noisy.data, FUN = pack_dwt)

#line one in first matrix of dwt_data contain wt transform of condition 1 in ind 1
plot(dwt_data[[1]][1,], c(wd(noisy.data[[1]][1, ])$D,wd(noisy.data[[1]][1, ])$C[length(wd(noisy.data[[1]][1, ])$C)]))

plot(dwt_data[[1]][2,], c(wd(noisy.data[[1]][2, ])$D,
                          wd(noisy.data[[1]][2, ])$C[length(wd(noisy.data[[1]][2, ])$C)]))







DW_tens <- rearrange( dwt_data, lev_res = lev_res, n_curve=3)


dim(DW_tens)

X <- G


#DW_tens[,1,] <- X[,1]%*%matrix(c(1,1,1), ncol=3) +rnorm(N*3, sd=0.1)
v1 <- rep( 1, N)
tt <- mv_reg  (DW_tens[,1,],cbind(v1,X[,1] ))# correspond to regression of covariate 1 on a n dimensional wavelet coefficient
tt
bmv_reg(DW_tens[,1,],cbind(v1,X[,1] ))


tt$Bhat/diag(tt$S)





### Fitting a mash per level of resolution

#idea generate two tensor, one coef one Shat
# then extract data to fit a mash

#need to parse by variable  (most inside loop) then by wave coef then by condition




Y <- DW_tens

marg_assoc <- cal_Bhat_Shat_tensor  (Y, X, v1)

marg_assoc
tens_marg <- marg_assoc

hist(marg_assoc $tens_Bhat[,,1]/marg_assoc $tens_Shat[,,1], nclass = 100)
image(marg_assoc$tens_Bhat[,,1]/marg_assoc $tens_Shat[,,1])


Bhat <- marg_assoc$tens_Bhat[,indx_lst[[3]],1]
Shat <- marg_assoc$tens_Shat[,indx_lst[[3]],1]
image(Bhat/Shat)

#fit mash

G_prior <- lapply(  1: (log2(dim(Y)[2])+1) , function(s) fit_mash_level(marg_assoc, s, indx_lst, data.driven = FALSE))

#or
G_prior <- init_prior_mvfsusie(tens_marg = tens_marg,
                               indx_list = indx_lst,
                               data.driven=FALSE)
class(G_prior)

#Now we need to compute the Bayes Factor to optimize the pi parameters
LBF <-  log_BF_tens  ( G_prior, tens_marg, indx_lst)



# the EM to maximize the marginal likelihood

res_EM <- EM_pi_mvfsusie(G_prior,
                         tens_marg,
                         indx_lst
)
EM_pi <- res_EM
#update the prior
Y <- DW_tens
mvfsusie_obj  <- init_mvfsusie_obj (L=4, G_prior, Y,X )
mvfsusie.obj <- mvfsusie_obj

G_prior <-  update_prior_weight_mvfsusie(G_prior,
                                         res_EM$tpi_k
)
class(G_prior)

alpha <-  cal_zeta_mvfsusie(  res_EM$lBF)
## Compute posterior
#How to reuse mash object
L=2
dim(mvfsusie_obj$fitted_wc[[1]])
dim(tens_marg$tens_Bhat)



post_tens <- get_post_tens (G_prior, tens_marg, indx_lst, all =TRUE)

str(post_tens)
#mvfsusie.obj <- update_pi(mvfsusie.obj, l=1, res_EM$tpi_k)
plot( update_pi( mvfsusie.obj =  mvfsusie.obj ,
                 l             = 1 ,
                 tpi           =  EM_pi$tpi_k)$est_pi[[1]][[5]],
      EM_pi$tpi_k[[5]]
)
plot( update_prior_weight_mvfsusie (get_G_prior(mvfsusie.obj)
                                    , EM_pi$tpi_k  )[[5]]$fitted_g$pi,
      EM_pi$tpi_k[[5]]
)


mvfsusie.obj <-  update_mvfsusie (mvfsusie.obj  = mvfsusie_obj ,
                                  l         = 1,
                                  EM_pi     = res_EM,
                                  tens_marg = tens_marg,
                                  indx_lst  = indx_lst,
                                  all=FALSE
)

plot( mvfsusie.obj$est_pi[[1]][[5]], res_EM$tpi_k[[5]])

plot( mvfsusie.obj$fitted_wc[[1]][,,1], post_tens$post_mean_tens[,,1])


mvfsusie.obj <-  update_mvfsusie (mvfsusie.obj  = mvfsusie_obj ,
                                  l         = 2,
                                  EM_pi     = res_EM,
                                  tens_marg = tens_marg,
                                  indx_lst  = indx_lst, all=TRUE
)


mvfsusie.obj$alpha





out <- mvfsusie(Y=noisy.data,
                X=X,
                L=2,maxit = 10)
out$alpha
out <- mvfsusie(Y=noisy.data,
                X=X,
                L=7,maxit = 30)
out$lBF
out$alpha
