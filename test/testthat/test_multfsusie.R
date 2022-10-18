
#### test with univariate and continuous data
library(testthat)
library(ashr)
library(wavethresh)
library(mixsqp)
library(susiF.alpha)
set.seed(1)
f1 <- simu_IBSS_per_level(lev_res=7, alpha=1, prop_decay =1.5)
# f1$true_coef effect is on wavelet coeff at scale 1 and 50 of coeff scale 2
plot(f1$sim_func, type="l", ylab="y")
set.seed(9)
f2 <- simu_IBSS_per_level(lev_res=6, alpha=0, prop_decay =1.5 )
# f2$true_coef effect is on wavelet coeff at scale 1 and 50 of coeff scale 2
plot(f2$sim_func, type="l", ylab="y")
N=500
P=10

nullweight=10
set.seed(23)
G = matrix(sample(c(0, 1,2), size=N*P, replace=T), nrow=N, ncol=P) #Genotype
beta0       <- 0
beta1       <- 1
pos1 <- 5
noisy.data  <- list()

low_wc=NULL
G = matrix(sample(c(0, 1,2), size=N*P, replace=T), nrow=N, ncol=P) #Genotype
beta1       <- 1
backfit=TRUE
greedy =TRUE
control_mixsqp=list(verbose=FALSE)

noisy.data  <- list()
 rsnr <- 6
for ( i in 1:N)
{
  f1_obs <- f1$sim_func
  f2_obs <- f2$sim_func
  Y0     <-  G[i,pos1]   +rnorm(1,sd=1)
  Y1     <- t(beta1*G[i,pos1]*f1_obs +  rnorm(length(f1$sim_func), sd=  (1/  rsnr ) *sd(f1$sim_func)))
  Y2     <- t(beta1*G[i,pos1]*f2_obs +  rnorm(length(f2$sim_func), sd=  (1/  rsnr ) *sd(f2$sim_func)))
  Y3     <- c(-1,1) *G[i, pos1]  + c(rnorm(2,sd=1))
  noisy.data [[i]] <- list(Y0,Y1,Y2,Y3)
}

X <- G

noisy.data[[1]]
noisy.data[[2]]
list_dfs  <- list()
for ( k in 1:length(noisy.data[[1]]))
{
  list_dfs [[k]]     <- do.call(rbind, lapply(noisy.data, `[[`, k))
}




plot(list_dfs[[2]][1,], type = "l", ylim = c(min(list_dfs[[2]]), max(list_dfs[[2]])))
for (i  in 2:N){

  lines(list_dfs[[2]][i,] , col=(G[i, pos1]+1))
}
plot(list_dfs[[3]][1,], type = "l", ylim = c(min(list_dfs[[3]]), max(list_dfs[[3]])))
for (i  in 2:N){

  lines(list_dfs[[3]][i,] , col=(G[i, pos1]+1))
}


type_mark <-  is.functional(list_dfs)

test_that("The inputs are of the following classes ",
          {
            expect_equal( type_mark$mark_type, c( "univariate", "functional" ,"functional" ,"univariate"))
            expect_equal( type_mark$dim_mark, c(  1 ,128 , 64  , 2))
            expect_equal( type_mark$ncond, 5)

          }
)


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

Y_data   <- list(Y_u =Y_u,
                 Y_f =Y_f)

X <- susiF.alpha:::colScale(X)
# centering input
Y_data <- multi_array_colScale(Y_data, scale=FALSE)
test_that("Y_u and Y_f should be not NULL ",
          {
            expect_equal( !is.null(Y_data$Y_u),TRUE)
            expect_equal( !is.null(Y_data$Y_f),TRUE)
          }
)




G_prior <- init_prior_multfsusie( Y=Y_data ,
                                 X=X,
                                 v1=v1,
                                 list_indx_lst=list_indx_lst,
                                 lowc_wc= low_wc,
                                 control_mixsqp=control_mixsqp,
                                 nullweight=  nullweight
)

test_that("G_prior object should have the following classes ",
          {
            expect_equal( class(G_prior),"multfsusie_prior")
            expect_equal( !is.null(G_prior$G_prior_u),TRUE)
            expect_equal( !is.null(G_prior$G_prior_f),TRUE)
            expect_equal( class(G_prior$G_prior_f[[1]]),"mixture_normal_per_scale")
            expect_equal( class(G_prior$G_prior_f[[2]]),"mixture_normal_per_scale")
          }
)
L=3
multfsusie.obj <- init_multfsusie_obj(L_max = 3, G_prior, Y_data,X , type_mark,
                                      greedy = greedy, backfit=backfit)
test_that("multfsusie internal prior to be equal to ",
          {
            multfsusie.obj <- init_multfsusie_obj(L_max = 3, G_prior, Y_data,X , type_mark,
                                                  greedy = greedy, backfit=backfit)

            expect_equal(get_G_prior (multfsusie.obj ),  G_prior)

          }
)
test_that("G_prior object should have the following classes ",
          {
            expect_equal( class(multfsusie.obj),"multfsusie")
            expect_equal( length(multfsusie.obj$fitted_wc),L)
            expect_equal( length(multfsusie.obj$fitted_wc2),L)
            expect_equal( length(multfsusie.obj$fitted_uni),L)
            expect_equal( length(multfsusie.obj$fitted_uni2),L)
            expect_equal( dim(multfsusie.obj$fitted_uni2[[1]] ),c(P, sum(type_mark$dim_mark[which(type_mark$mark_type=="univariate")] ) ))
            expect_equal( dim(multfsusie.obj$fitted_uni[[1]] ),c(P, sum(type_mark$dim_mark[which(type_mark$mark_type=="univariate")] ) ))
            expect_equal( length(multfsusie.obj$fitted_wc[[1]]), length(which(type_mark$mark_type=="functional")))
            expect_equal( length(multfsusie.obj$fitted_wc2[[1]]), length(which(type_mark$mark_type=="functional")))
          }
)

update_Y <-Y_data
effect_estimate   <- cal_Bhat_Shat_multfsusie(update_Y,X,v1)




test_that("The estimated effect should be", {
  expect_equal(effect_estimate$res_f[[1]]$Bhat[pos1,-128], f1$true_coef, tol= 2*(1/sqrt(N)))# removing C coefficient
  expect_equal(effect_estimate$res_f[[2]]$Bhat[pos1,-64], f2$true_coef, tol= 2*(1/sqrt(N)))# removing C coefficient
  expect_equal(effect_estimate$res_uni $Bhat[pos1, ], c(1,-1,1), tol= 2*(1/sqrt(N)))

})



test_that("G_prior is correctly updated",
          {
            tpi               <- get_pi(multfsusie.obj,1)
            tpi$est_pi_f[[1]][[1]][1] <- 10
            G_prior <- update_prior(G_prior, tpi= tpi)
            expect_equal(get_pi_G_prior(G_prior)$est_pi_f[[1]][[1]][1] , 10)
          }
  )
tpi               <- get_pi(multfsusie.obj,1)
G_prior <- update_prior(G_prior, tpi= tpi) #allow EM to start close to previous solution (to double check)
init_pi0_w= 1
control_mixsqp =  list(
                        eps = 1e-6,
                        numiter.em = 40,
                        verbose = FALSE
                      )

EM_out  <- EM_pi_multsusie(G_prior  = G_prior,
                           effect_estimate= effect_estimate,
                           list_indx_lst =  list_indx_lst,
                           init_pi0_w= init_pi0_w,
                           control_mixsqp =  control_mixsqp
)

zeta <- cal_zeta(EM_out$lBF)
test_that("The highest assignation should be equal to", {
  zeta <- cal_zeta(EM_out$lBF)
  expect_equal(which.max(zeta), pos1
  )
})

tpi <- EM_out$tpi_k
test_that("The output class should be", {
expect_equal(class(tpi$est_pi_f[[1]]),"pi_mixture_normal_per_scale")
expect_equal(class(tpi$est_pi_f[[2]]),"pi_mixture_normal_per_scale")
expect_equal(class(tpi$est_pi_u[[1]]),"pi_mixture_normal")
expect_equal(class(tpi$est_pi_u[[2]]),"pi_mixture_normal")
expect_equal(class(tpi$est_pi_u[[3]]),"pi_mixture_normal")
})

L_mat <- L_mixsq_multsusie (G_prior, effect_estimate, list_indx_lst)


test_that("The output class should be", {
  expect_equal(class(L_mat),"lik_multfsusie")
  expect_equal(!is.null(L_mat$L_mat_u),TRUE)
  expect_equal(!is.null(L_mat$L_mat_f),TRUE)
})

est_pi_f <- lapply(1:length(L_mat$L_mat_f) ,
                   function(k) m_step(L_mat$L_mat_f[[k]],
                                      zeta,
                                      list_indx_lst[[k]],
                                      init_pi0_w= init_pi0_w,
                                      control_mixsqp = control_mixsqp
                   )
)


test_that("The highest assignation should be equal to", {
  expect_equal(tpi$est_pi_u[[1]][1], 0,
               tolerance = 0.01)
  expect_equal(tpi$est_pi_u[[2]][1], 0,
               tolerance = 0.01)
  expect_equal(tpi$est_pi_u[[3]][1], 0,
               tolerance = 0.01)
  expect_lt( get_pi0(tpi = tpi$est_pi_f[[1]])[1 ], c(0.1  ) )
  expect_lt( get_pi0(tpi = tpi$est_pi_f[[1]])[2], c( 0.6 ) )
  expect_lt( get_pi0(tpi = tpi$est_pi_f[[2]])[1 ], c(0.1  ) )
  expect_lt( get_pi0(tpi = tpi$est_pi_f[[2]])[2], c( 0.6 ) )
})



multfsusie.obj <- update_multfsusie(multfsusie.obj  = multfsusie.obj ,
                                    l               = 1,
                                    EM_pi           = EM_out,
                                    effect_estimate = effect_estimate,
                                    list_indx_lst   = list_indx_lst)
