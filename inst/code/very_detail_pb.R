rm(list=ls())
devtools::load_all(".")

source("D:/Document/Serieux/Travail/Package/mvf.susie.alpha/R/computational_routine.R")


load("D:/Document/Serieux/Travail/Package/mvf.susie.alpha/inst/pb_input.RData")
X <- pb_input[[3]]$X
Y <- pb_input[[3]]$Y
true_pos  <- pb_input[[3]]$true_pos

source("D:/Document/Serieux/Travail/Package/mvf.susie.alpha/R/operation_on_multfsusie_obj.R")




L=11
data.format = "list_df"#"ind_mark"
verbose=TRUE
maxit = 100
tol = 1e-3
cov_lev = 0.95
min.purity=0.5
L_start=3
#all = FALSE
filter.cs =TRUE
init_pi0_w=1
nullweight =10
control_mixsqp =  list(
  eps = 1e-6,
  numiter.em = 40,
  verbose = FALSE
)
thresh_lowcount=0
cal_obj=FALSE
greedy=TRUE
backfit=TRUE
parallel=FALSE
max_SNP_EM=1000
gridmult=sqrt(2)
max_step_EM=1
cor_small=FALSE

thresh_lowcount=0

#### formating -----




if(verbose){
  print("Data transform")
}

ind_analysis <- which_notNA_pos(Y)
h <- 1
list_wdfs <- list()
list_indx_lst  <-  list()
if( !is.null(Y$Y_f)){
  outing_grid <- list()


  pos <- list()
  for (i in 1:length(Y$Y_f))
  {
    pos[[i]] <- 1:ncol(Y$Y_f[[i]])
  }

  for (i in 1:length(Y$Y_f)){
    if ( !(length(pos[[i]])==ncol(Y$Y_f[[i]]))) #miss matching positions and number of observations
    {
      stop(paste("Error: number of position provided different from the number of column of Y$Y_f, entry",i))
    }
  }



  for ( k in 1:length(Y$Y_f))
  {
    map_data <-  susiF.alpha::remap_data(Y=Y$Y_f[[k]],
                                         pos=pos[[k]],
                                         verbose=verbose)
    outing_grid[[k]] <- map_data$outing_grid

    temp               <- DWT2( map_data$Y)
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

type_mark <- is.functional ( Y=Y_data,
                             data.format=data.format)







#### centering and scaling covariate ----
X <- susiF.alpha:::colScale(X)

# centering input
Y_data <- multi_array_colScale(Y_data, scale=FALSE)
#
if(verbose){
  print("Data transform done")
}

### Cleaning ------

#### discarding  null/low variance    ------

threshs <- create_null_thresh(type_mark = type_mark)

low_trait <- check_low_count  (Y_data, thresh_lowcount=threshs,ind_analysis = ind_analysis  )

v1 <- rep( 1, nrow(X))
if(verbose){
  print("Initializing prior")
}
### Work here ------

if( cor_small){
  df <- list()

  if( !is.null(ind_analysis$idx_u)){
    df$Y_u <- lengths(ind_analysis$idx_u)-1
  }else{
    df$Y_u <- NULL
  }
  if( !is.null(ind_analysis$idx_f)){
    df$Y_f <- lengths(ind_analysis$idx_f)-1
  }else{
    df$Y_f <- NULL
  }

}else{
  df= NULL
}

temp  <- init_prior_multfsusie(Y              = Y_data ,
                               X              = X,
                               v1             = v1,
                               list_indx_lst  = list_indx_lst,
                               low_trait      = low_trait,
                               control_mixsqp = control_mixsqp,
                               nullweight     = nullweight,
                               ind_analysis   = ind_analysis,
                               parallel       = parallel,
                               max_SNP_EM     = max_SNP_EM,
                               gridmult       = gridmult,
                               max_step_EM    = max_step_EM
)

G_prior          <- temp$G_prior
effect_estimate  <- temp$res
init             <- TRUE


multfsusie.obj <- init_multfsusie_obj( L_max=L,
                                       G_prior=G_prior,
                                       Y=Y_data,
                                       X=X,
                                       type_mark=type_mark,
                                       L_start=L_start,
                                       greedy=greedy,
                                       backfit=backfit,
                                       ind_analysis   = ind_analysis)


set.seed(1)

for (i in 1:length(multfsusie.obj$fitted_wc)){


  for ( j in 1: length(multfsusie.obj$fitted_wc[[i]]) ){
    multfsusie.obj$fitted_wc[[i]][[j]] <- matrix(rnorm(n=prod(dim(multfsusie.obj$fitted_wc[[i]][[j]]))),
                                                 nrow= nrow(multfsusie.obj$fitted_wc[[i]][[j]]))
  }
}


check <- 3*tol

update_Y    <-  Y_data

#### start while -----

iter <- 1   ####ITER 1------


 l=1
  ####l 1------
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

  multfsusie.obj <- update_multfsusie(multfsusie.obj  = multfsusie.obj ,
                                      l               = l,
                                      EM_pi           = EM_out,
                                      effect_estimate = effect_estimate,
                                      list_indx_lst   = list_indx_lst,
                                      low_trait       = low_trait )




  multfsusie.obj$alpha[[1]] <- 0*  multfsusie.obj$alpha[[1]]
  multfsusie.obj$alpha[[1]][207] <- 1


####l 2------

l=2
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

  multfsusie.obj <- update_multfsusie(multfsusie.obj  = multfsusie.obj ,
                                      l               = l,
                                      EM_pi           = EM_out,
                                      effect_estimate = effect_estimate,
                                      list_indx_lst   = list_indx_lst,
                                      low_trait       = low_trait )







  multfsusie.obj$alpha

  ####l 3------

  l=3
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

  multfsusie.obj <- update_multfsusie(multfsusie.obj  = multfsusie.obj ,
                                      l               = l,
                                      EM_pi           = EM_out,
                                      effect_estimate = effect_estimate,
                                      list_indx_lst   = list_indx_lst,
                                      low_trait       = low_trait )







  multfsusie.obj$alpha

  temp <- wavethresh::wd(rep(0, 32))


  temp$D     <- multfsusie.obj$fitted_wc[[3]][[1]][388,-32]
  temp$C[length(temp$C)]     <- multfsusie.obj$fitted_wc[[3]][[1]][388, 32]
  plot(wavethresh::wr(temp),type="l")



  temp$D     <- multfsusie.obj$fitted_wc[[3]][[1]][410,-32]
  temp$C[length(temp$C)]     <- multfsusie.obj$fitted_wc[[3]][[1]][410, 32]
  plot(wavethresh::wr(temp),type="l")











multfsusie.obj <- update_alpha_hist(multfsusie.obj)

multfsusie.obj <- update_cal_cs(multfsusie.obj,
                                cov_lev=cov_lev)
multfsusie.obj$cs
#caca <- multfsusie.obj
#multfsusie.obj  <- caca

dummy.cs <-  which_dummy_cs(multfsusie.obj,
                            min.purity = min.purity,
                            X=X)


dummy.cs

multfsusie.obj <- greedy_backfit (multfsusie.obj,
                                  verbose        = verbose,
                                  cov_lev        = cov_lev,
                                  X              = X,
                                  min.purity     = min.purity
)

sigma2    <- estimate_residual_variance(multfsusie.obj,
                                        Y=Y_data,
                                        X=X,
                                        ind_analysis = ind_analysis)
multfsusie.obj <- update_residual_variance(multfsusie.obj,
                                           sigma2 = sigma2 )

multfsusie.obj <- test_stop_cond(multfsusie.obj      = multfsusie.obj,
                                 check               = check,
                                 cal_obj             = cal_obj,
                                 Y                   = Y_data,
                                 X                   = X,
                                 list_indx_lst       = list_indx_lst,
                                 ind_analysis        = ind_analysis)
check <- multfsusie.obj$check



iter <- iter+1



#end while
