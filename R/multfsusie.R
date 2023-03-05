#' @title Sum of multiple Single Functions
#'
#' @description Implementation of the multfSuSiE method
#'
#' @details Implementation of the multfSuSiE method
#'
#' @param Y list of observed time series. Length of N in which every element
#' contains a xi (number of condition) by 2^S matrix. The matrix corresponds to the
#' individuals multivariate time series
#' @param X matrix of size n by p contains the covariates
#'@param pos  list of sampling position for Y$Y_f entries, if not provided the alogirthm will consider that the column are venly spaced
#'if only one of the entry (say entry 2 out of 3 entries) has some non evenly spaced data then you can defined pos as follow,
#'pos =list(pos1=NULL,pos2=pos_uneven, pos3=NULL), where pos_uneven  is a vector containning the sampling position (assumed to be increasing)
#')
#' @param L the number of effect to fit (if not specified set to =2)
#'
#' @param pos vector of length J, corresponding to position/time pf
#' the observed column in each entery of Y$Y_f, if missing suppose that the observation
#' are evenly spaced
#' @param L_start number of effect initialized at the start of the algorithm
#'@param data.format character specify hw the input data is organised,
#' "ind_mark" the input is a list in which each element is a list of individual mark measurement.
#'  "list_df", corresponds to the case where the input is a list of  of data frames
#'   in which element from univariate trait are stored in Y$Y_u, one column corresponds to a univariate trait
#'    (can be set to NULL if no univariate trait considered) and functional trait are stored in the sub list Y$Y_f
#'    where each element of the sub list  Y$Y_f is a n by T data frame (T being the number of observation points)
#'    (can be NULL if no functional trait considered)
#' @param control_mixsqp list of parameter for mixsqp function see  mixsqp package
#'
#' @param verbose If \code{verbose = TRUE}, the algorithm's progress,
#' and a summary of the optimization settings, are printed to the
#' console.
#'
#' @param thresh_lowcount numeric, use to check the wavelet coefficients have
#'  problematic distribution (very low dispersion even after standardization).
#'  Basically check if the median of the absolute value of the distribution of
#'   a wavelet coefficient is below this threshold, if yes the algorithm discard
#'   this wavelet coefficient (setting its estimate effect to 0 and estimate sd to 1).
#'   Set to 0 by default. Can be useful when analyzing sparse data from sequence
#'    based assay or small samples.
#' @param greedy logical, if true allow greedy search for extra effect
#'  (up to L specify by the user). Set as TRUE by default
#' @param backfit logical, if true allow discarding effect via backfitting.
#'  Set as true by default as TRUE. We advise to keep it as TRUE
#'
#' @param tol A small, non-negative number specifying the convergence
#' tolerance for the IBSS fitting procedure. The fitting procedure
#' will halt when the difference in the variational lower bound, or
#' \dQuote{ELBO} (the objective function to be maximized), is less
#' than \code{tol}.
#'
#' @param maxit Maximum number of IBSS iterations to perform.
#'
#' @param cov_lev numeric between 0 and 1, corresponding to the
#' expected level of coverage of the cs if not specified set to 0.95
#' @param  cal_obj logical if set as true compute ELBO for convergence monitoring
#' @param min.purity minimum purity for estimated credible sets
#' @param filter.cs logical, if TRUE filter the credible set (removing low purity cs and cs with estimated prior equal to 0)
#'@param thresh_lowcount list  of numeric (check example), use to check the
#'   wavelet coefficients/univariate trait have with
#'  problematic distribution (very low dispersion even after standardization).
#'  Basically check if the median of the absolute value of the distribution of
#'   a wavelet coefficient/univariate trait is below a given threshold, if yes the algorithm discard
#'   this wavelet coefficient (setting its estimate effect to 0 and estimate sd to 1).
#'   Set to 0 by default. Can be useful when analyzing sparse data from sequence
#'    based assay or small samples.
#' @param nullweight numeric value for penalizing likelihood at point mass 0 (should be between 0 and 1)
#' (usefull in small sample size)
#' @param init_pi0_w starting value of weight on null compoenent in mixsqp
#'  (between 0 and 1)
#' @export
#' @examples
#'
#'
#'library(mvf.susie.alpha)
#'set.seed(1)
#'
#'N <- 100 #Sample size
#'P= 100 # number of SNP
#'L <- sample(1:10, size=1) #Number of effect
#'print(L)
#'list_lev_res <- list(5,6) # two functional phenotypes ,
#'#one of length 2^5, and one of length 2^6)
#'n_univ <- 3 #3 univariate phenotypes
#'eff <-  list()
#'for(l in 1:L){ #Simulate the mult-trait effect
#'  eff[[l]] <-   simu_effect_multfsusie (list_lev_res=list_lev_res,
#'                                        n_univ=n_univ, output_level = 2)
#'}
#'
#'
#'Y_f1 <-  matrix(rnorm((2^list_lev_res[[1]])*N ,sd=1), nrow = N)
#'Y_f2 <-  matrix(rnorm((2^list_lev_res[[2]])*N ,sd=1), nrow = N)
#'
#'Y_u <- matrix(rnorm((n_univ)*N ,sd=1), nrow = N)
#'
#'
#'G = matrix(sample(c(0, 1,2), size=N*P, replace=TRUE), nrow=N, ncol=P) #Genotype
#'
#'
#'true_pos <- sample( 1:ncol(G), L)# actually causal column/SNP
#'
#'for ( i in 1:N){
#'  for ( l in 1:L){
#'
#'    Y_f1[i,]<- Y_f1[i,]+eff[[l]]$func_effect[[1]]$sim_func*G[i,true_pos[[l]]]
#'    Y_f2[i,]<- Y_f2[i,]+eff[[l]]$func_effect[[2]]$sim_func*G[i,true_pos[[l]]]
#'    Y_u[i,]<- Y_u[i,]+ eff[[l]]$univ_effect*G[i,true_pos[[l]]]
#'  }
#'}
#'
#'Y_f <- list()
#'Y_f[[1]] <- Y_f1
#'Y_f[[2]] <- Y_f2
#'Y <- list( Y_f = Y_f, Y_u=Y_u) # preparing data ,
#'#current ouput type expect list of two which element named
#'#Y_f for functional trait and Y_u for univariate trait
#'
#'
#'pos = list(pos1= 1: ncol(Y$Y_f[[1]]),
#'           pos2= 1: ncol(Y$Y_f[[2]])) # if you signal is sample between 1 and 64
#'
#'m1 <- multfsusie(Y=Y,
#'                 X=G,
#'                 pos=pos,
#'                 L=11 ,
#'                 data.format="list_df",
#'                 L_start=11 ,
#'                 nullweight=10,
#'                 maxit=10)
#'m1$cs# credible sets
#'
#'
#'##thresholding some trait
#'
#'#create object for trhesholding some trait for a user specified value thresholding
#'
#'threshs <- threshold_set_up( thresh_u= rep(1e-3,3), thresh_f = c(1e-3, 1e-3))
#'
#'m1 <- multfsusie(Y=Y,
#'                 X=G,
#'                 L=11 ,
#'                 data.format="list_df",
#'                 L_start=11,
#'                 thresh_lowcount=threshs,
#'                 maxit=10)

multfsusie <- function(Y ,X,L=2, pos = NULL,
                       data.format = "list_df",#"ind_mark",
                       verbose=TRUE,
                       maxit = 100,
                       tol = 1e-3,
                       cov_lev = 0.95,
                       min.purity=0.5,
                       L_start=3,
                       #all = FALSE,
                       filter.cs =TRUE,
                       init_pi0_w=1,
                       nullweight ,
                       control_mixsqp =  list(
                                              eps = 1e-6,
                                              numiter.em = 40,
                                              verbose = FALSE
                                             ),
                       thresh_lowcount,
                       cal_obj=FALSE,
                       greedy=TRUE,
                       backfit=TRUE
                      )


{
  pt <- proc.time()
  if(missing(nullweight )){
    nullweight <- 10#/sqrt(nrow(X))
  }

  if(L_start >L)
  {
    L_start <- L
  }
#Formatting the data ----
####ind mark type ----
  if(data.format=="ind_mark")  {
  list_dfs  <- list()
  for ( k in 1:length(Y[[1]]))
  {
    list_dfs [[k]]     <- do.call(rbind, lapply(Y, `[[`, k))
  }

  type_mark <-  is.functional(list_dfs)

  list_wdfs <- list()
  list_indx_lst  <-  list()
  if( "functional" %!in% type_mark$mark_type)
  {
    Y_f <- NULL
  }else{

    if(missing(pos)){
      pos <- list()
      for ( k in which(type_mark$mark_type=="functional"))
      {
        pos[[k]] <- 1:ncol(list_dfs[[k]])
      }
    }
    for (i in 1:length(which(type_mark$mark_type=="functional"))){
      if ( !(length(pos[[i]])==ncol(list_dfs[[k]]))) #miss matching positions and number of observations
      {
        stop(paste("Error: number of position provided different from the number of column of Y$Y_f, entry",i))
      }
    }




    outing_grid <- list()
    h <- 1
    for ( k in which(type_mark$mark_type=="functional"))
    {

      map_data <-  susiF.alpha::remap_data(Y=list_dfs[[k]],
                                           pos=pos[[k]],
                                           verbose=verbose)
      outing_grid[[k]] <- map_data$outing_grid

      temp               <- DWT2(map_data$Y)
      list_wdfs[[h]]     <- cbind( temp$D,temp$C)
      list_indx_lst[[h]] <- susiF.alpha::gen_wavelet_indx( log2(ncol(  list_wdfs[[h]]) ))
      h <- h+1
      rm(map_data)
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
}####list_df ----
  if(data.format=="list_df"){

    h <- 1
    list_wdfs <- list()
    list_indx_lst  <-  list()
    if( !is.null(Y$Y_f)){
      outing_grid <- list()

      if(missing(pos)){
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

  }

  #### centering and scaling covariate ----
  X <- susiF.alpha:::colScale(X)

  # centering input
  Y_data <- multi_array_colScale(Y_data, scale=FALSE)
  #

  ### Cleaning ------

  #### discarding  null/low variance    ------
  if( missing(thresh_lowcount)){
    threshs <- create_null_thresh(type_mark = type_mark)
  }
  low_trait <- check_low_count  (Y_data, thresh_lowcount=threshs  )

  v1 <- rep( 1, nrow(X))

  temp  <- init_prior_multfsusie(Y              = Y_data ,
                                 X              = X,
                                 v1             = v1,
                                 list_indx_lst  = list_indx_lst,
                                 low_trait      = low_trait,
                                 control_mixsqp = control_mixsqp,
                                 nullweight     = nullweight
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
                                         backfit=backfit)


  # numerical value to check breaking condition of while
  check <- 3*tol

  update_Y    <-  Y_data


  if( L==1)
  {

    effect_estimate   <- cal_Bhat_Shat_multfsusie(update_Y,X,v1,low_trait=low_trait)
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
                                                 list_indx_lst  = indx_lst
                                  )
    )

    sigma2    <- estimate_residual_variance(multfsusie.obj,Y=Y_data,X)
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
                                                       low_trait = low_trait)
         tpi               <- get_pi(multfsusie.obj,l)
         G_prior <- update_prior(G_prior, tpi= tpi ) #allow EM to start close to previous solution (to double check)

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
                                            low_trait       = low_trait)


      }#end for l in 1:L  -----

      multfsusie.obj <- greedy_backfit (multfsusie.obj,
                                        verbose        = verbose,
                                        cov_lev        = cov_lev,
                                        X              = X,
                                        min.purity     = min.purity
                                        )

      sigma2         <- estimate_residual_variance.multfsusie(multfsusie.obj,Y=Y_data,X)
      multfsusie.obj <- update_residual_variance(multfsusie.obj, sigma2 = sigma2 )

      multfsusie.obj <- test_stop_cond(multfsusie.obj = multfsusie.obj,
                                  check               = check,
                                  cal_obj             = cal_obj,
                                  Y                   = Y_data,
                                  X                   = X,
                                  list_indx_lst       = list_indx_lst)
      check <- multfsusie.obj$check



      iter <- iter+1



    }#end while
  }


  #preparing output
   multfsusie.obj <- out_prep(multfsusie.obj  = multfsusie.obj,
                                   Y          = Y_data,
                                   X          = X,
                              list_indx_lst   = list_indx_lst,
                                   filter.cs  = filter.cs,
                                 outing_grid  = outing_grid
    )
   multfsusie.obj$runtime <- proc.time()-pt
  return(multfsusie.obj)

}
