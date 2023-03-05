library(susiF.alpha)
library(mvf.susie.alpha)
N <- 100 #Sample size
P= 100 # number of SNP
L <- sample(1:10, size=1) #Number of effect
print(L)
list_lev_res <- list(5,6) # two functional phenotypes ,
#one of length 2^5, and one of length 2^6)
n_univ <- 3 #3 univariate phenotypes
eff <-  list()
for(l in 1:L){ #Simulate the mult-trait effect
  eff[[l]] <-   simu_effect_multfsusie (list_lev_res=list_lev_res,
                                        n_univ=n_univ, output_level = 2)
}


Y_f1 <-  matrix(rnorm((2^list_lev_res[[1]])*N ,sd=1), nrow = N)
Y_f2 <-  matrix(rnorm((2^list_lev_res[[2]])*N ,sd=1), nrow = N)

Y_u <- matrix(rnorm((n_univ)*N ,sd=1), nrow = N)


G = matrix(sample(c(0, 1,2), size=N*P, replace=TRUE), nrow=N, ncol=P) #Genotype


true_pos <- sample( 1:ncol(G), L)# actually causal column/SNP

for ( i in 1:N){
  for ( l in 1:L){

    Y_f1[i,]<- Y_f1[i,]+eff[[l]]$func_effect[[1]]$sim_func*G[i,true_pos[[l]]]
    Y_f2[i,]<- Y_f2[i,]+eff[[l]]$func_effect[[2]]$sim_func*G[i,true_pos[[l]]]
    Y_u[i,]<- Y_u[i,]+ eff[[l]]$univ_effect*G[i,true_pos[[l]]]
  }
}


#### ADDing NAs------
Y_f <- list()
Y_f[[1]] <- Y_f1
Y_f[[1]][1:10,]<-NA

Y_f[[2]] <- Y_f2


Y_u[50:60,1]<-NA

Y_u[55:70,2]<-NA
Y <- list( Y_f = Y_f, Y_u=Y_u) # preparing data ,
#current ouput type expect list of two which element named
#Y_f for functional trait and Y_u for univariate trait

Y_f1 -Y_f[[1]]
Y_f[[1]][is.na(Y_f[[1]])]<-0
Y_f[[1]]

greedy=TRUE
backfit=TRUE
L_start=4
init_pi0_w=0.9
nullweight=10
control_mixsqp=list(verbose=FALSE)
#### START ALGO------
which_notNA_pos <-  function( Y){

  if( !is.null(Y$Y_f)){
    temp_f <-  lapply( 1:length(Y$Y_f), function(d)
                                          which(complete.cases(Y$Y_f[[d]]))
      )

  }else{
    temp_f <- NULL
  }
  if( !is.null(Y$Y_u)){
    temp_u <-    lapply( 1:ncol(Y$Y_u), function(d)
                                               which(complete.cases(Y$Y_u[,d]))

                       )

    if(length(temp_u)==0){
      temp_u <- NULL
    }


  }else{
    temp_u <- NULL
  }
  out <- list( idx_f =temp_f,
               idx_u  =  temp_u)
  return( out)
}

ind_analysis <- which_notNA_pos(Y)
pos = list(pos1= 1: ncol(Y$Y_f[[1]]),
           pos2= 1: ncol(Y$Y_f[[2]])) # if you signal is sample between 1 and 64
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
                             data.format="list_df")

X <- G

X <- susiF.alpha:::colScale(X)

# centering input
Y_data <- multi_array_colScale(Y_data, scale=FALSE)
### good up to here ----




#### discarding  null/low variance    ------
#if( missing(thresh_lowcount)){
  threshs <- create_null_thresh(type_mark = type_mark)
#}
low_trait <- check_low_count  (Y_data, thresh_lowcount=threshs,ind_analysis = ind_analysis  )

v1 <- rep( 1, nrow(X))

temp  <- init_prior_multfsusie(Y              = Y_data ,
                               X              = X,
                               v1             = v1,
                               list_indx_lst  = list_indx_lst,
                               low_trait      = low_trait,
                               control_mixsqp = control_mixsqp,
                               nullweight     = nullweight,
                               ind_analysis   = ind_analysis
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


check <- 3*tol

update_Y    <-  Y_data


if( L==1)
{

  effect_estimate   <- cal_Bhat_Shat_multfsusie(update_Y,X,v1,
                                                low_trait=low_trait,
                                                ind_analysis   = ind_analysis)
  tpi               <- get_pi(multfsusie.obj,1)
  G_prior           <- update_prior(G_prior, tpi= tpi) #allow EM to start close to previous solution (to double check)


  EM_out  <- EM_pi_multsusie(G_prior         = G_prior,
                             effect_estimate = effect_estimate,
                             list_indx_lst   = list_indx_lst,
                             init_pi0_w      = init_pi0_w,
                             control_mixsqp  = control_mixsqp,
                             nullweight      = nullweight,
                             low_trait       = low_trait
  )



  multfsusie.obj <- update_multfsusie(multfsusie.obj  = multfsusie.obj ,
                                      l               = 1,
                                      EM_pi           = EM_out,
                                      effect_estimate = effect_estimate,
                                      list_indx_lst   = list_indx_lst,
                                      low_trait       = low_trait )


  multfsusie.obj <- update_ELBO(multfsusie.obj,
                                get_objective( multfsusie.obj = multfsusie.obj,
                                               Y         = Y_data ,
                                               X         = X,
                                               list_indx_lst  = indx_lst,
                                               ind_analysis = ind_analysis
                                )
  )

  sigma2    <- estimate_residual_variance(multfsusie.obj,
                                          Y=Y_data,
                                          X=X,
                                          ind_analysis = ind_analysis)
  multfsusie.obj <- update_residual_variance(multfsusie.obj, sigma2 = sigma2 )

}else{
  ##### Start While -----
  iter <- 1
  while(check >tol & (h/L) <maxit)
  {

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
                                                      low_trait=low_trait,
                                                      ind_analysis   = ind_analysis)
        tpi               <- get_pi(multfsusie.obj,1)
        G_prior           <- update_prior(G_prior, tpi= tpi) #allow EM to start close to previous solution (to double check)


        EM_out  <- EM_pi_multsusie(G_prior         = G_prior,
                                   effect_estimate = effect_estimate,
                                   list_indx_lst   = list_indx_lst,
                                   init_pi0_w      = init_pi0_w,
                                   control_mixsqp  = control_mixsqp,
                                   nullweight      = nullweight,
                                   low_trait       = low_trait
        )


      }

      multfsusie.obj <- update_multfsusie(multfsusie.obj  = multfsusie.obj ,
                                          l               = l,
                                          EM_pi           = EM_out,
                                          effect_estimate = effect_estimate,
                                          list_indx_lst   = list_indx_lst,
                                          low_trait       = low_trait )


    }#end for l in 1:L  -----

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



  }#end while
}


