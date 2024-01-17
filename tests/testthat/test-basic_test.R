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
f3 <- simu_IBSS_per_level(lev_res=5, alpha=0, prop_decay =1.5 )

# f2$true_coef effect is on wavelet coeff at scale 1 and 50 of coeff scale 2
plot(f2$sim_func, type="l", ylab="y")
N=100
P=20

nullweight= 0
set.seed(23)
G = matrix(sample(c(0, 1,2), size=N*P, replace=T), nrow=N, ncol=P) #Genotype
beta0       <- 0
beta1       <- 1
pos1 <- 5
pos2 <- 10
pos3 <- 1
noisy.data  <- list()

low_wc=NULL
G = matrix(sample(c(0, 1,2), size=N*P, replace=T), nrow=N, ncol=P) #Genotype
beta1          <- 1
backfit        = TRUE
greedy         = TRUE
verbose        = TRUE
df             = NULL
parallel       = FALSE
control_mixsqp = list(verbose=FALSE)
cov_lev        = 0.95
min.purity     = 0.5
filter.number = 10
family = "DaubLeAsymm"
pos=NULL
prior = "mixture_normal"

noisy.data  <- list()
rsnr <- 6
for ( i in 1:N)
{
  f1_obs <- f1$sim_func
  f2_obs <- f2$sim_func
  f3_obs <- f3$sim_func

  Y0     <-  G[i,pos1]   +rnorm(1,sd=1)
  Y1     <- t(beta1*G[i,pos1]*f1_obs +  beta1*G[i,pos3]*f3_obs +rnorm(length(f1$sim_func), sd=  (1/  rsnr ) *sd(f1$sim_func)))
  Y2     <- t(beta1*G[i,pos2]*f2_obs +  rnorm(length(f2$sim_func), sd=  (1/  rsnr ) *sd(f2$sim_func)))
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


Y <- list( Y_u= cbind( list_dfs[[1]], list_dfs[[4]]),
           Y_f= list_dfs[ c(2,3 ) ])

plot(Y$Y_f[[1]][1,], type = "l", ylim = c(min(Y$Y_f[[1]]), max(Y$Y_f[[1]])))
for (i  in 2:N){

  lines(Y$Y_f[[1]][i,] , col=(G[i, pos1]+1))
}
plot(Y$Y_f[[2]][1,], type = "l", ylim = c(min(Y$Y_f[[2]]), max(Y$Y_f[[2]])))
for (i  in 2:N){

  lines(Y$Y_f[[2]][i,] , col=(G[i, pos1]+1))
}


type_mark <-  is.functional(Y )

test_that("The inputs are of the following classes ",
          {
            expect_equal( type_mark$mark_type, c(  "functional" ,"functional" ,"univariate"))
            expect_equal( type_mark$dim_mark, c(   128 , 64  , 3))
            expect_equal( type_mark$ncond, 5)

          }
)


h <- 1
list_wdfs <- list()
list_indx_lst  <-  list()
if( !is.null(Y$Y_f)){
  outing_grid <- list()

  if( is.null(pos)){
    pos <- list()
    for (i in 1:length(Y$Y_f))
    {
      pos[[i]] <- 1:ncol(Y$Y_f[[i]])
    }
  }
  for (i in 1:length(Y$Y_f)){
    if ( !(length(pos[[i]])==ncol(Y$Y_f[[i]]))) #miss matching positions and number of observations
    {
      stop(paste("Error: number of position provided different from the number of column of Y$Y_f, entry",i))
    }
  }


  interpolated_Y <- Y


  for ( k in 1:length(Y$Y_f))
  {
    map_data <-  susiF.alpha::remap_data(Y=Y$Y_f[[k]],
                                         pos=pos[[k]],
                                         verbose=verbose)
    outing_grid[[k]] <- map_data$outing_grid
    interpolated_Y$Y_f[[k]] <-  map_data$Y



    temp               <- DWT2( map_data$Y,
                                filter.number = filter.number,
                                family        = family)
    list_wdfs[[h]]     <- cbind( temp$D,temp$C)
    list_indx_lst[[h]] <- susiF.alpha:::gen_wavelet_indx( log2(ncol(  list_wdfs[[h]]) ))
    h <- h+1
    rm(map_data)
  }
  Y_f <- list_wdfs
  v1  <- nrow( Y_f [[1]])
}else{
  Y_f <- NULL
  v1  <- nrow( Y_u)
}


Y_data   <- list(Y_u =Y$Y_u,
                 Y_f =Y_f)



X_old <- X
Y_old <-Y_data
X <- susiF.alpha:::colScale(X)

# centering input
Y_data <- multi_array_colScale(Y_data, scale=FALSE)
test_that("Correct variance estimates  for X component",{
  expect_equal( attr(X,"scaled:center"),apply(X_old,2, mean))
  expect_equal( attr(X,"scaled:scale"),apply(X_old,2,sd))
  expect_equal( attr(X,"scaled:scale"),apply(X_old,2,sd))
  expect_equal( attr(Y_data$Y_u,"scaled:center"),apply(Y_old$Y_u,2, mean))
  expect_equal( attr(Y_data$Y_f[[1]],"scaled:center"),apply(Y_old$Y_f[[1]],2, mean))
  expect_equal( attr(Y_data$Y_f[[2]],"scaled:center"),apply(Y_old$Y_f[[2]],2, mean))

}
)


ind_analysis <- which_notNA_pos(Y_data)

test_that("The number of individuals for each mark should be",
          {
            Y_data2 <-Y_data
            idx_na <- sample (1:nrow(Y_data2$Y_u), size=1)
            Y_data2$Y_u[idx_na,2] <-NA

            expect_equal(which_notNA_pos(Y_data2)$idx_u[[2]] ,(1:nrow(Y_data2$Y_u))[-idx_na])
            idx_na <- sample (1:nrow(Y_data2$Y_u), size=1)
            Y_data2$Y_f[[1]][idx_na] <-NA
            expect_equal(which_notNA_pos(Y_data2)$idx_f[[1]] ,(1:nrow(Y_data2$Y_f[[1]]))[-idx_na])

          }


          )

threshs <- create_null_thresh(type_mark = type_mark)
low_trait <- check_low_count  (Y_data, thresh_lowcount=threshs  )



test_that("Checking if the null traits are captured",
          {
            Y_data2 <-Y_data
            Y_data2$Y_u[,2]<-0* Y_data2$Y_u[,2]

            Y_data2$Y_f[[2]][,1:20]<-0* Y_data2$Y_f[[2]][,1:20]
            threshs <- create_null_thresh(type_mark = type_mark)
            low_trait <- check_low_count  (Y_data2, thresh_lowcount=threshs  )
            expect_equal(low_trait$low_u,2)
            expect_equal(low_trait$low_wc[[2]],c(1:20))
            expect_equal(low_trait$low_wc[[1]],NULL)
          }


)


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
                                  low_trait= low_trait,
                                  control_mixsqp=control_mixsqp,
                                  nullweight=  nullweight,
                                  max_SNP_EM = 100
)$G_prior

test_that("G_prior object should have the following classes ",
          {
            expect_equal( class(G_prior),"multfsusie_prior")
            expect_equal( !is.null(G_prior$G_prior_u),TRUE)
            expect_equal( !is.null(G_prior$G_prior_f),TRUE)
            expect_equal( class(G_prior$G_prior_f[[1]]),"mixture_normal")
            expect_equal( class(G_prior$G_prior_f[[2]]),"mixture_normal")
          }
)
L=3
multfsusie.obj <- init_multfsusie_obj(L_max          = 10,
                                      G_prior        = G_prior,
                                      Y              = Y_data,
                                      X              = X,
                                      type_mark      = type_mark,
                                      L_start        =   3,
                                      greedy         = greedy,
                                      backfit        = backfit,
                                      ind_analysis   = ind_analysis)
test_that("multfsusie internal prior to be equal to ",
          {
            multfsusie.obj <- init_multfsusie_obj(L_max          = 10,
                                                  G_prior        = G_prior,
                                                  Y              = Y_data,
                                                  X              = X ,
                                                  type_mark      = type_mark,
                                                  L_start        =   3,
                                                  greedy         = greedy,
                                                  backfit        = backfit,
                                                  ind_analysis   = ind_analysis)
            expect_equal(get_G_prior (multfsusie.obj ),  G_prior)

          }
)
test_that("G_prior object should have the following classes ",
          {
            expect_equal( class(multfsusie.obj),"multfsusie")
            expect_equal( length(multfsusie.obj$fitted_wc),L)
            expect_equal( length(multfsusie.obj$fitted_wc2),L)
            expect_equal( length(multfsusie.obj$fitted_u),L)
            expect_equal( length(multfsusie.obj$fitted_u2),L)
            expect_equal( dim(multfsusie.obj$fitted_u2[[1]] ),c(P, sum(type_mark$dim_mark[which(type_mark$mark_type=="univariate")] ) ))
            expect_equal( dim(multfsusie.obj$fitted_u[[1]] ),c(P, sum(type_mark$dim_mark[which(type_mark$mark_type=="univariate")] ) ))
            expect_equal( length(multfsusie.obj$fitted_wc[[1]]), length(which(type_mark$mark_type=="functional")))
            expect_equal( length(multfsusie.obj$fitted_wc2[[1]]), length(which(type_mark$mark_type=="functional")))
          }
)

update_Y <-Y_data
effect_estimate   <- cal_Bhat_Shat_multfsusie(update_Y,X,v1)




test_that("The estimated effect should be", {
  expect_equal(effect_estimate$res_f[[1]]$Bhat[pos1,-128]/attr(X,"scaled:scale")[pos1], f1$true_coef,
               tolerance = 2*(1/sqrt(N)))# removing C coefficient
  expect_equal(effect_estimate$res_f[[2]]$Bhat[pos1,-64]/attr(X,"scaled:scale")[pos1], f2$true_coef,
               tolerance= 2*(1/sqrt(N)))# removing C coefficient
  expect_equal(effect_estimate$res_u $Bhat[pos1, ]/attr(X,"scaled:scale")[pos1], c(1,-1,1),
               tolerance= 2*(1/sqrt(N)))

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

EM_out  <- EM_pi_multsusie(G_prior         = G_prior,
                           effect_estimate = effect_estimate,
                           list_indx_lst   = list_indx_lst,
                           init_pi0_w      = init_pi0_w,
                           control_mixsqp  = control_mixsqp,
                           nullweight      = nullweight,
                           low_trait       = low_trait
)

zeta <- cal_zeta(EM_out$lBF)
test_that("The highest assignation should be equal to", {
  zeta <- cal_zeta(EM_out$lBF)
  expect_equal(which.max(zeta), pos3
  )
})

tpi <- EM_out$tpi_k
test_that("The output class should be", {
  expect_equal(class(tpi$est_pi_f[[1]]),"pi_mixture_normal")
  expect_equal(class(tpi$est_pi_f[[2]]),"pi_mixture_normal")
  expect_equal(class(tpi$est_pi_u[[1]]),"pi_mixture_normal")
  expect_equal(class(tpi$est_pi_u[[2]]),"pi_mixture_normal")
  expect_equal(class(tpi$est_pi_u[[3]]),"pi_mixture_normal")
})

L_mat <- L_mixsq_multsusie (G_prior, effect_estimate, list_indx_lst, idx=1:ncol(X))


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
                                      control_mixsqp = control_mixsqp,
                                      nullweight = nullweight
                   )
)


test_that("The highest assignation should be equal to", {
  expect_equal(tpi$est_pi_u[[1]][1], 0,
               tolerance = 0.01)
  expect_equal(tpi$est_pi_u[[2]][1], 0,
               tolerance = 0.01)
  expect_equal(tpi$est_pi_u[[3]][1], 0,
               tolerance = 0.01)
  expect_lt( susiF.alpha::get_pi0(tpi = tpi$est_pi_f[[1]]),  c(0.73  ) )
  expect_lt( susiF.alpha::get_pi0(tpi = tpi$est_pi_f[[2]]) , c( 999998 ) )
})

threshs <- create_null_thresh(type_mark = type_mark)
low_trait <- check_low_count  (Y_data, thresh_lowcount=threshs  )


test_that("check greedy backfit",{

  iter=1
  init=TRUE
  for( l in 1:multfsusie.obj$L)
  {

    update_Y <- cal_partial_resid(multfsusie.obj = multfsusie.obj,
                                  l              = (l-1)  ,
                                  X              = X,
                                  Y              = Y_data,
                                  list_indx_lst  = list_indx_lst
    )



    if(verbose){
      print(paste("Fitting effect ", l,", iter" ,  iter ))
    }
    if(init){#recycle operation used to fit the prior

      EM_out <- susiF.alpha:::gen_EM_out (tpi_k= get_pi_G_prior(G_prior),
                                          lBF  = log_BF(G_prior,
                                                        effect_estimate,
                                                        list_indx_lst,
                                                        low_trait = low_trait)
      )
      class(EM_out) <- c("EM_pi_multfsusie","list")
      init <- FALSE
    }else{

      effect_estimate   <- cal_Bhat_Shat_multfsusie(update_Y,X,v1,
                                                    low_trait      = low_trait,
                                                    ind_analysis   = ind_analysis,
                                                    parallel       = parallel
      )
      tpi               <- get_pi(multfsusie.obj,1)
      G_prior           <- update_prior(G_prior, tpi= tpi) #allow EM to start close to previous solution (to double check)


      EM_out  <- EM_pi_multsusie(G_prior         = G_prior,
                                 effect_estimate = effect_estimate,
                                 list_indx_lst   = list_indx_lst,
                                 init_pi0_w      = init_pi0_w,
                                 control_mixsqp  = control_mixsqp,
                                 nullweight      = nullweight,
                                 low_trait       = low_trait,
                                 df              = df
      )


    }

  multfsusie.obj <- update_multfsusie(multfsusie.obj  = multfsusie.obj ,###TODO:: SLOW
                                      l               = l,
                                      EM_pi           = EM_out,
                                      effect_estimate = effect_estimate,
                                      list_indx_lst   = list_indx_lst,
                                      low_trait       = low_trait )


  }#end for l in 1:L  -----

  multfsusie.obj$lfsr_wc
  multfsusie.obj$lfsr_u
  multfsusie.obj <- greedy_backfit (multfsusie.obj,
                                    verbose        = verbose,
                                    cov_lev        = cov_lev,
                                    X              = X,
                                    min.purity     = min.purity
  )
  iter =  iter+1
  expect_equal(multfsusie.obj$L ,  (7+3))

  expect_equal(length(multfsusie.obj$fitted_wc),  (7+3))
  expect_equal(length(multfsusie.obj$fitted_u),  (7+3))
  expect_equal( multfsusie.obj$cs[[1]], 1)
  expect_equal( multfsusie.obj$cs[[2]], 5)
  expect_equal( multfsusie.obj$cs[[3]], 10)


  for( l in 1:multfsusie.obj$L)
  {

    update_Y <- cal_partial_resid(multfsusie.obj = multfsusie.obj,
                                  l              = (l-1)  ,
                                  X              = X,
                                  Y              = Y_data,
                                  list_indx_lst  = list_indx_lst
    )



    if(verbose){
      print(paste("Fitting effect ", l,", iter" ,  iter ))
    }
    if(init){#recycle operation used to fit the prior

      EM_out <- susiF.alpha:::gen_EM_out (tpi_k= get_pi_G_prior(G_prior),
                                          lBF  = log_BF(G_prior,
                                                        effect_estimate,
                                                        list_indx_lst,
                                                        low_trait = low_trait)
      )
      class(EM_out) <- c("EM_pi_multfsusie","list")
      init <- FALSE
    }else{

      effect_estimate   <- cal_Bhat_Shat_multfsusie(update_Y,X,v1,
                                                    low_trait      = low_trait,
                                                    ind_analysis   = ind_analysis,
                                                    parallel       = parallel
      )
      tpi               <- get_pi(multfsusie.obj,1)
      G_prior           <- update_prior(G_prior, tpi= tpi) #allow EM to start close to previous solution (to double check)


      EM_out  <- EM_pi_multsusie(G_prior         = G_prior,
                                 effect_estimate = effect_estimate,
                                 list_indx_lst   = list_indx_lst,
                                 init_pi0_w      = init_pi0_w,
                                 control_mixsqp  = control_mixsqp,
                                 nullweight      = nullweight,
                                 low_trait       = low_trait,
                                 df              = df
      )


    }

    multfsusie.obj <- update_multfsusie(multfsusie.obj  = multfsusie.obj ,###TODO:: SLOW
                                        l               = l,
                                        EM_pi           = EM_out,
                                        effect_estimate = effect_estimate,
                                        list_indx_lst   = list_indx_lst,
                                        low_trait       = low_trait )


  }#end for l in 1:L  -----

  multfsusie.obj$lfsr_wc
  multfsusie.obj$lfsr_u
  multfsusie.obj <- greedy_backfit (multfsusie.obj,
                                    verbose        = verbose,
                                    cov_lev        = cov_lev,
                                    X              = X,
                                    min.purity     = min.purity
  )

  multfsusie.obj$cs
})



##### la -----

test_that("The performance  of multfsusie on  this example should be",{
  library(susiF.alpha)
  library(mvf.susie.alpha)
  set.seed(1)
  N=100
  P=50
  G = matrix(sample(c(0, 1,2), size=N*P, replace=T), nrow=N, ncol=P) #Genotype
  beta1       <- 1
  beta2       <- 1
  L <- 4#actual number of effect
  lf <-  list()
  for(l in 1:L){
    lf[[l]] <- simu_IBSS_per_level(lev_res=5)$sim_func #functional effect for effect l
  }

  tt <- sample(0:4,1)
  if( length(which(apply(G,2,var)==0))>0){
    G <- G[,-which(apply(G,2,var)==0)]
  }
  # G <- matrix( rnorm(100*300), nrow = 100)
  true_pos <- sample( 1:ncol(G), L)

  Y <- matrix(rnorm((2^5)*100 ,sd=4), nrow = 100)
  for ( i in 1:100){
    for ( l in 1:L){
      Y[i,] <- Y[i,]+ lf[[l]]*G[i,true_pos[[l]]]
    }
  }
  Y_f <- list()
  Y_f[[1]] <- Y

  Y <- list( Y_f = Y_f, Y_u=NULL)

  m1 <- multfsusie(Y=Y,
                   X=G,
                   L=11 ,
                   L_start=11 ,
                   nullweight=10,
                   cal_obj =FALSE,
                   maxit=10, verbose=FALSE)

  expect_lt(  sqrt(mean( (m1$fitted_func[[2]][[1]]-lf[[3]])^2)),0.7)
  expect_lt(  sqrt(mean( (m1$fitted_func[[1]][[1]]-lf[[1]])^2)),0.7)
  expect_lt(  sqrt(mean( (m1$fitted_func[[3]][[1]]-lf[[2]])^2)),0.7)
  expect_equal( length(which(true_pos%in% do.call(c, m1$cs))) , length(true_pos))
}
)


##### lfsr -----

G_prior <- get_G_prior(multfsusie.obj)

cal_clfsr(G_prior         = get_G_prior(multfsusie.obj),
          effect_estimate = effect_estimate,
          list_indx_lst   = list_indx_lst)





