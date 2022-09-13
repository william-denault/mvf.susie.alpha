#rm(list = ls())
library(ashr)
library(mashr)
library(abind)
library(wavethresh)
library(mvtnorm)
library(mixsqp)
library(susiF.alpha)
library(mvf.susie.alpha)
set.seed(13)
L=3

min_levres <- 4

#gene expression effect of SNP 1 and SNP2
f01  <-  1
f02  <- -0


#effect of SNP 1 and SNP2 on mark 1
lev_res1 <- 5
temp_func <-  simu_IBSS_per_level(lev_res1)
f11 <-  temp_func$sim_func
temp_func <-  simu_IBSS_per_level(lev_res1)
f12 <-  0*temp_func$sim_func
#effect of SNP 1 and SNP2 on mark 2
lev_res2 <- 7
temp_func <-  simu_IBSS_per_level(lev_res2)
f21 <-  temp_func$sim_func
temp_func <-  simu_IBSS_per_level(lev_res2)
f22 <- 0* temp_func$sim_func

#effect of SNP 1 and SNP2 on 2 CpGs

f31  <-  c(1,-1)
f32  <-  0*c(-1,1)





indx_lst1 <- susiF.alpha::gen_wavelet_indx(lev_res = lev_res1)
indx_lst2 <- susiF.alpha::gen_wavelet_indx(lev_res = lev_res2)

N = 100

#Number of covariates

P = 10

#Choosing which variable will have an effect
pos1 <- 1
pos2 <- 2



G = matrix(sample(c(0, 1,2), size=N*P, replace=T), nrow=N, ncol=P) #Genotype
beta1       <- 1



noisy.data  <- list()
for ( i in 1:N)
{

  Y0               <- f01*G[i,pos1]+f02*G[i,pos2] +rnorm(1,sd=0.08) #gene expression
  Y1               <- t( beta1*G[i,pos1]*f11+ beta1*G[i,pos2]*f12+rnorm(2^lev_res1 ,sd=0.08))#mark 1
  Y2               <- t( beta1*G[i,pos1]*f21+ beta1*G[i,pos2]*f22 +rnorm(2^lev_res2 ,sd=0.08))#mark 2
  Y3               <-  f31*G[i, pos1]+f32*G[i, pos2] + c(rnorm(2,sd=0.08))
  noisy.data [[i]] <-  list(Y0,Y1,Y2,Y3)

}

noisy.data[[1]]
noisy.data[[2]]
list_dfs  <- list()
for ( k in 1:length(noisy.data[[1]]))
{
  list_dfs [[k]]     <- do.call(rbind, lapply(noisy.data, `[[`, k))
}

type_mark <-  is.functional(list_dfs)
type_mark

list_wdfs <- list()
list_indx_lst  <-  list()
if( "functional" %!in% type_mark$mark_type)
{
  Y_f <- NULL
}else{
  h <- 1
  for ( k in which(type_mark$mark_type=="functional"))
  {
    temp               <- DWT2(list_dfs[[k]])
    list_wdfs[[h]]     <- cbind( temp$D,temp$C)
    list_indx_lst[[h]] <- gen_wavelet_indx( log2(ncol(  list_wdfs[[h]]) ))
    h <- h+1
  }
  Y_f <- list_wdfs
  v1  <- nrow( Y_f [[1]])
}
if("univariate" %!in% type_mark$mark_type)
{
  Y_u <- NULL
}else{
  Y_u <- do.call( cbind, list_dfs [ which(type_mark$mark_type=="univariate") ])
  v1  <- nrow(Y_u)
}

X   <- G
Y   <- list(Y_u =Y_u,
            Y_f =Y_f)

G_prior <- init_prior_multfsusie(Y,
                                 X,
                                 v1,
                                 list_indx_lst
                                 )


multfsusie.obj <- init_multfsusie_obj(L, G_prior, Y,X , type_mark)
class(multfsusie.obj)
update_Y    <-  Y



  l=1

  # numerical value to check breaking condition of while
  check <- 1
  h     <- 0

  tt   <- cal_Bhat_Shat_multfsusie(update_Y,X,v1)
  tpi  <- get_pi(multfsusie.obj ,1)
  tpi$est_pi_u[[1]][1] <-16
  G_prior <- update_prior(G_prior, tpi= tpi ) #allow EM to start close to previous solution (to double check)

  class(G_prior$G_prior_f[[1]])
  class(G_prior$G_prior_u[[1]])
  G_prior$G_prior_u[[1]]$fitted_g$pi[[1]]==16 #good lord it works

  effect_estimate= tt


  EM_out  <- EM_pi_multsusie(G_prior  = G_prior,
                             effect_estimate= tt,
                             list_indx_lst =  list_indx_lst,
                             control_mixsqp =  list(
                               eps = 1e-3,
                               numiter.em = 40,
                               verbose = FALSE
                             )
  )

  EM_pi <-  EM_out
  EM_pi$lBF



  multfsusie.obj <- update_multfsusie(multfsusie.obj   = multfsusie.obj ,
                                      l               = l,
                                      EM_pi           = EM_out,
                                      effect_estimate = effect_estimate,
                                      list_indx_lst   = list_indx_lst)
  multfsusie.obj$fitted_uni[[l]]
  multfsusie.obj$fitted_wc[[l]]

  update_Y <- cal_partial_resid(multfsusie.obj = multfsusie.obj,
                                l              = l ,
                                X              = X,
                                Y              = Y,
                                list_indx_lst  = list_indx_lst
  )
  multfsusie.obj$lBF
  multfsusie.obj$alpha

  plot( update_Y$Y_u[,1], X[,1])
  plot( Y$Y_u[,1], X[,1])

  l=2

  # numerical value to check breaking condition of while
  check <- 1
  h     <- 0

  tt   <- cal_Bhat_Shat_multfsusie(update_Y,X,v1)
  tpi  <- get_pi(multfsusie.obj ,1)
  tpi$est_pi_u[[1]][1] <-16
  G_prior <- update_prior(G_prior, tpi= tpi ) #allow EM to start close to previous solution (to double check)

  class(G_prior$G_prior_f[[1]])
  class(G_prior$G_prior_u[[1]])
  G_prior$G_prior_u[[1]]$fitted_g$pi[[1]]==16 #good lord it works

  effect_estimate= tt


  EM_out  <- EM_pi_multsusie(G_prior  = G_prior,
                             effect_estimate= tt,
                             list_indx_lst =  list_indx_lst,
                             control_mixsqp =  list(
                               eps = 1e-3,
                               numiter.em = 40,
                               verbose = FALSE
                             )
  )

  EM_pi <-  EM_out
  EM_pi$lBF



  multfsusie.obj <- update_multfsusie(multfsusie.obj   = multfsusie.obj ,
                                      l               = l,
                                      EM_pi           = EM_out,
                                      effect_estimate = effect_estimate,
                                      list_indx_lst   = list_indx_lst)
  multfsusie.obj$lBF
  multfsusie.obj$alpha

  multfsusie.obj$fitted_uni[[l]]
  multfsusie.obj$fitted_wc[[l]]

  update_Y <- cal_partial_resid(multfsusie.obj = multfsusie.obj,
                                l              = l ,
                                X              = X,
                                Y              = Y,
                                list_indx_lst  = list_indx_lst
  )




  l=3

  # numerical value to check breaking condition of while
  check <- 1
  h     <- 0

  tt   <- cal_Bhat_Shat_multfsusie(update_Y,X,v1)
  tpi  <- get_pi(multfsusie.obj ,1)
  tpi$est_pi_u[[1]][1] <-16
  G_prior <- update_prior(G_prior, tpi= tpi ) #allow EM to start close to previous solution (to double check)

  class(G_prior$G_prior_f[[1]])
  class(G_prior$G_prior_u[[1]])
  G_prior$G_prior_u[[1]]$fitted_g$pi[[1]]==16 #good lord it works

  effect_estimate= tt


  EM_out  <- EM_pi_multsusie(G_prior  = G_prior,
                             effect_estimate= tt,
                             list_indx_lst =  list_indx_lst,
                             control_mixsqp =  list(
                               eps = 1e-3,
                               numiter.em = 40,
                               verbose = FALSE
                             )
  )

  EM_pi <-  EM_out
  EM_pi$lBF



  multfsusie.obj <- update_multfsusie(multfsusie.obj   = multfsusie.obj ,
                                      l               = l,
                                      EM_pi           = EM_out,
                                      effect_estimate = effect_estimate,
                                      list_indx_lst   = list_indx_lst)
  multfsusie.obj$lBF
  multfsusie.obj$alpha

  multfsusie.obj$fitted_uni[[l]]
  multfsusie.obj$fitted_wc[[l]]

  update_Y <- cal_partial_resid(multfsusie.obj = multfsusie.obj,
                                l              = l ,
                                X              = X,
                                Y              = Y,
                                list_indx_lst  = list_indx_lst
  )

  multfsusie.obj$lBF
  multfsusie.obj$alpha
  plot(multfsusie.obj$alpha[[1]], col="red")
  points(multfsusie.obj$alpha[[2]], col="blue")
  points(multfsusie.obj$alpha[[3]], col="green")



multfsusie.obj$est_pi[[1]]

multfsusie.obj$est_pi[[2]]

multfsusie.obj$est_pi[[3]]


library(mvf.susie.alpha)
out <- multfsusie(Y=noisy.data,
                X=G,
                L=6,maxit = 20)
out$alpha



plot( out$ELBO)

