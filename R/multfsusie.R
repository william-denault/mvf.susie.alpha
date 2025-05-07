#' @title Multi-modal fSuSiE
#' @description Implementation of the Multi-modal fSuSiE method.
#'
#' @details Implementation of the multfSuSiE method.
#'
#' @param Y A list of data frames. Univariate traits are stored in Y$Y_u, with one column per univariate trait (can be set to NULL if no univariate traits are considered). Functional traits are stored in the sublist Y$Y_f, where each element of Y$Y_f is an n by T data frame (T being the number of observation points) (can be NULL if no functional trait is considered).
#'
#' @param X A matrix of size n by p containing the covariates.
#'
#' @param pos A list of sampling positions for Y$Y_f entries.
#' If not provided, the algorithm will assume the columns are evenly spaced.
#'  If only one entry (e.g., entry 2 out of 3) has non-evenly spaced data,
#'   you can define pos as follows: pos = list(pos1=NULL, pos2=pos_uneven, pos3=NULL),
#'    where pos_uneven is a vector containing the sampling positions (assumed to be increasing).
#'
#' @param L The number of effects to fit (default is 2).
#'
#' @param post_processing Character, use "TI" for translation-invariant wavelet estimates,
#'  "HMM" for hidden Markov model (useful for estimating non-zero regions),
#'  or "none" for simple wavelet estimate (not recommended).
#'
#' @param prior Specifies the prior used in functional trait analysis.
#'  The two available choices are "mixture_normal_per_scale" and "mixture_normal" (default).
#'   Using "mixture_normal" is up to 40\% faster but may lead to slight power loss.
#'
#' @param L_start Number of effects initialized at the start of the algorithm.
#'
#' @param control_mixsqp A list of parameters for the mixsqp function (see mixsqp package).
#'
#' @param verbose If TRUE, the algorithm's progress and a summary of the optimization
#' settings are printed to the console.
#'
#' @param thresh_lowcount Numeric, used to check if wavelet coefficients have
#' problematic distribution (very low dispersion even after standardization).
#' If the median of the absolute value of the distribution of a wavelet coefficient
#'  is below this threshold, the algorithm discards this wavelet coefficient
#'   (setting its estimated effect to 0 and its estimated sd to 1).
#'   Set to 0 by default. Useful when analyzing sparse data from sequence-based
#'    assays or small samples.
#'
#' @param greedy Logical, if TRUE allows greedy search for extra effects
#'  (up to L specified by the user). Set to TRUE by default.
#'
#' @param backfit Logical, if TRUE allows discarding effects via backfitting.
#'  Set to TRUE by default. It is advised to keep this as TRUE.
#'
#' @param tol A small, non-negative number specifying the convergence tolerance
#' for the IBSS fitting procedure. The fitting procedure will halt when
#'  the difference in the variational lower bound, or ELBO
#'   (the objective function to be maximized), is less than tol.
#'
#' @param maxit Maximum number of IBSS iterations to perform.
#'
#' @param cov_lev Numeric between 0 and 1, corresponding to the expected
#' level of coverage of the credible sets. Default is 0.95.
#'
#' @param cal_obj Logical, if TRUE computes ELBO for convergence monitoring.
#'
#' @param min_purity Minimum purity for estimated credible sets.
#'
#' @param filter_cs Logical, if TRUE filters the credible sets
#' (removing low purity credible sets and those with estimated prior equal to 0).
#'
#' @param nullweight Numeric value for penalizing likelihood at point mass 0
#' (should be between 0 and 1). Useful in small sample sizes.
#'
#' @param min_purity minimum purity for estimated credible sets
#'
#' @param filter_cs logical, if TRUE filters the credible set (removing low-purity)
#' cs and cs with estimated prior equal to 0). Set as TRUE by default.
#'
#' @param gridmult numeric used to control the number of components used in the mixture prior (see ashr package
#'  for more details). From the ash function:  multiplier by which the default grid values for mixed differ from one another.
#'   (Smaller values produce finer grids.). Increasing this value may reduce computational time.
#'
#' @param max_scale numeric, define the maximum of wavelet coefficients used in the analysis (2^max_scale).
#'        Set 10 true by default.
#'
#' @param max_step_EM max_step_EM
#'
#' @param max_SNP_EM maximum number of SNP used for learning the prior. By default, set to 1000. Reducing this may help reduce
#' the computational time. We advise to keep it at least larger than 50
#'
#' @param cor_small logical set to FALSE by default. If TRUE used the Bayes factor from Valen E Johnson JRSSB 2005 instead of Wakefield approximation for Gen Epi 2009
#' The Bayes factor from Valen E Johnson JRSSB 2005 tends to have better coverage in small sample sizes. We advise using this parameter if n<50
#'
#' @param filter.number see documentation of wd from wavethresh package
#'
#' @param family see documentation of wd from wavethresh package
#'
#' @param init_pi0_w Starting value of weight on null
#' component in mixsqp (between 0 and 1).
#'
#' @param e Threshold value to avoid computing posterior probabilities
#'  that have low alpha values. Set to 0 to compute the entire posterior.
#'  Default value is 0.001.
#'
#' @export
#' @examples
#' library(mvf.susie.alpha)
#' set.seed(1)
#'
#' N <- 100 # Sample size
#' P <- 100 # Number of SNPs
#' L <- 2 # Number of effects
#' list_lev_res <- list(5, 6)
#' # Two functional phenotypes, one of length 2^5 and one of length 2^6
#' n_univ <- 3 # 3 univariate phenotypes
#' eff <- list()
#' for (l in 1:L) { # Simulate the multi-trait effect
#'   eff[[l]] <- simu_effect_multfsusie(list_lev_res=list_lev_res,
#'                                      n_univ=n_univ, output_level=2)
#' }
#'
#' Y_f1 <- matrix(rnorm((2^list_lev_res[[1]]) * N, sd=1), nrow=N)
#' Y_f2 <- matrix(rnorm((2^list_lev_res[[2]]) * N, sd=1), nrow=N)
#' Y_u <- matrix(rnorm(n_univ * N, sd=1), nrow=N)
#' # Genotype matrix
#' G <- matrix(sample(c(0, 1, 2), size=N*P, replace=TRUE), nrow=N, ncol=P)
#'
#' true_pos <- sample(1:ncol(G), L) # Actually causal columns/SNPs
#'
#' for (i in 1:N) {
#'   for (l in 1:L) {
#'     Y_f1[i,] <- Y_f1[i,] + eff[[l]]$func_effect[[1]]$sim_func * G[i, true_pos[[l]]]
#'     Y_f2[i,] <- Y_f2[i,] + eff[[l]]$func_effect[[2]]$sim_func * G[i, true_pos[[l]]]
#'     Y_u[i,]  <- Y_u[i,]  + eff[[l]]$univ_effect * G[i, true_pos[[l]]]
#'   }
#' }
#'
#' Y_f <- list(Y_f1, Y_f2)
#' Y <- list(Y_f=Y_f, Y_u=Y_u) # Preparing data
#'
#' pos <- list(pos1=1:ncol(Y$Y_f[[1]]), pos2=1:ncol(Y$Y_f[[2]]))
#' # If your signal is sampled between 1 and 64
#'
#' m1 <- multfsusie(Y=Y, X=G, pos=pos, L=3)
#' print(m1$cs) # Credible sets
#' print(true_pos)
#'
#' par(mfrow=c(2,1))
#' plot(m1$fitted_func[[1]][[1]], type="l", col="green",
#' main="Estimated function for the first marker", ylab="y")
#' lines(eff[[2]]$func_effect[[1]]$sim_func)
#' lines(m1$cred_band[[1]][[1]][1,], lty=2, col="darkgreen")
#' lines(m1$cred_band[[1]][[1]][2,], lty=2, col="darkgreen")
#'
#' plot(m1$fitted_func[[1]][[2]], type="l", col="green",
#' main="Estimated function for the second marker", ylab="y")
#' lines(eff[[2]]$func_effect[[2]]$sim_func)
#' lines(m1$cred_band[[1]][[2]][1,], lty=2, col="darkgreen")
#' lines(m1$cred_band[[1]][[2]][2,], lty=2, col="darkgreen")
#'
#' # A bit slower but useful for properly estimating the support of the effect
#' m1 <- multfsusie(Y=Y, X=G, pos=pos, L=3, post_processing = "HMM")
#'
#' plot(m1$fitted_func[[1]][[1]], type="l", col="green",
#'  main="Estimated function for the first marker", ylab="y")
#' lines(eff[[2]]$func_effect[[1]]$sim_func)
#' abline(h=0)
#' lines(m1$lfsr[[1]]$est_lfsr_functional[[1]], lty=2, col="darkgreen")
#' plot(m1$fitted_func[[1]][[2]], type="l", col="green",
#' main="Estimated function for the first marker", ylab="y")
#' lines(eff[[2]]$func_effect[[2]]$sim_func)
#' abline(h=0)
#' lines(m1$lfsr[[1]]$est_lfsr_functional[[2]], lty=2, col="darkgreen")
multfsusie <- function(Y, X, L = 2,
                       pos = NULL,
                       prior = "mixture_normal",
                       post_processing=c("smash","TI","HMM","none"),
                       verbose = TRUE,
                       maxit = 100,
                       tol = 1e-3,
                       cov_lev = 0.95,
                       min_purity = 0.5,
                       L_start = 3,
                       filter_cs = TRUE,
                       init_pi0_w = 1,
                       nullweight = 10,
                       control_mixsqp = list(
                         eps = 1e-6,
                         numiter.em = 40,
                         verbose = FALSE
                       ),
                       thresh_lowcount,
                       cal_obj = FALSE,
                       greedy = TRUE,
                       backfit = TRUE,
                       max_SNP_EM = 100,
                       gridmult = sqrt(2),
                       max_scale = 10,
                       max_step_EM = 1,
                       cor_small = FALSE,
                       filter.number = 10,
                       family = "DaubLeAsymm",
                       e = 0.001,
                       tol_null_prior=0.001)


{
  pt <- proc.time()

  prior           <- match.arg(prior)
  post_processing <- match.arg( post_processing)
  if(L_start >L)
  {
    L_start <- L
  }
  if(verbose){
    print("Starting initialization")
  }

  ind_analysis <- which_notNA_pos(Y)
#remove column of X constant in some sub cases
  tidx <- check_cst_X_sub_case(X,ind_analysis)


  if( length(tidx)>0){
    warning(paste("Some of the columns of X are constants, we removed" ,length(tidx), "columns"))
    X <- X[,-tidx]
  }
#Formatting the data ----

    if(verbose){
      print("Data transform")
    }


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
        map_data <-  fsusieR::remap_data(Y=Y$Y_f[[k]],
                                             pos=pos[[k]],
                                             verbose=verbose)
        outing_grid[[k]] <- map_data$outing_grid
        interpolated_Y$Y_f[[k]] <-  map_data$Y



        temp               <- DWT2( map_data$Y,
                                    filter.number = filter.number,
                                    family        = family)
        list_wdfs[[h]]     <- cbind( temp$D,temp$C)
        list_indx_lst[[h]] <- fsusieR::gen_wavelet_indx( log2(ncol(  list_wdfs[[h]]) ))
        h <- h+1
        rm(map_data)
      }
      Y_f <- list_wdfs
      v1  <- nrow( Y_f [[1]])
    }else{
      Y_f <- NULL
      v1  <- nrow( Y$Y_u)
      outing_grid <- NULL
    }


    Y_data   <- list(Y_u =Y$Y_u,
                     Y_f =Y_f)

    type_mark <- is.functional ( Y=Y_data  )



  #### centering and scaling covariate ----
  X <- fsusieR::colScale(X)

  # centering input
  Y_data <- multi_array_colScale(Y_data, scale=FALSE)
  #
  if(verbose){
    print("Data transform done")
  }

  ### Cleaning ------

  #### discarding  null/low variance    ------
  if( missing(thresh_lowcount)){
    threshs <- create_null_thresh(type_mark = type_mark)
  }else{
    threshs <- thresh_lowcount
  }
  low_trait <- check_low_count  (Y_data,
                                 thresh_lowcount = threshs,
                                 ind_analysis    = ind_analysis
                                 )

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
                                 prior          = prior,
                                 list_indx_lst  = list_indx_lst,
                                 low_trait      = low_trait,
                                 control_mixsqp = control_mixsqp,
                                 nullweight     = nullweight,
                                 ind_analysis   = ind_analysis,
                                 max_SNP_EM     = max_SNP_EM,
                                 gridmult       = gridmult,
                                 max_step_EM    = max_step_EM,
                                 tol_null_prior = tol_null_prior
  )

  G_prior          <- temp$G_prior
  effect_estimate  <- temp$res
  init             <- TRUE


  multfsusie.obj <- init_multfsusie_obj( L_max         = L,
                                         G_prior       = G_prior,
                                         Y             = Y_data,
                                         X             = X,
                                         type_mark     = type_mark,
                                         L_start       = L_start,
                                         greedy        = greedy,
                                         backfit       = backfit,
                                         ind_analysis  = ind_analysis)


  check <- 3*tol

 # browser()
  update_Y    <-  Y_data

  if(verbose){
    print("Initialization done")
  }
  if( L==1)
  {

    effect_estimate   <- cal_Bhat_Shat_multfsusie(update_Y,X,v1,
                                                  low_trait=low_trait,
                                                  ind_analysis   = ind_analysis )
    tpi               <- get_pi(multfsusie.obj,1)
    G_prior           <- update_prior(G_prior, tpi= tpi) #allow EM to start close to previous solution (to double check)


    EM_out  <- EM_pi_multsusie(G_prior         = G_prior,
                               effect_estimate = effect_estimate,
                               list_indx_lst   = list_indx_lst,
                               init_pi0_w      = init_pi0_w,
                               control_mixsqp  = control_mixsqp,
                               nullweight      = nullweight,
                               low_trait       = low_trait,
                               max_SNP_EM      = max_SNP_EM,
                               max_step        = max_step_EM,
                               df              = df
    )



    multfsusie.obj <- update_multfsusie(multfsusie.obj  = multfsusie.obj ,
                                        l               = 1,
                                        EM_pi           = EM_out,
                                        effect_estimate = effect_estimate,
                                        list_indx_lst   = list_indx_lst,
                                        low_trait       = low_trait ,
                                        e               = e)


    multfsusie.obj <- update_ELBO(multfsusie.obj,
                                  get_objective( multfsusie.obj = multfsusie.obj,
                                                 Y         = Y_data ,
                                                 X         = X,
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
    while(check >tol & iter <=maxit)
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

          EM_out <- fsusieR:::gen_EM_out (tpi_k= get_pi_G_prior(G_prior),
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
                                                        ind_analysis   = ind_analysis
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
                                            low_trait       = low_trait,
                                            e               = e)


      }#end for l in 1:L  -----

      multfsusie.obj <- greedy_backfit (multfsusie.obj,
                                        verbose        = verbose,
                                        cov_lev        = cov_lev,
                                        X              = X,
                                        min_purity     = min_purity
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
  }#end else in if(L==1)



#browser()
  #preparing output
   multfsusie.obj <- out_prep(multfsusie.obj  = multfsusie.obj,
                              Y               = Y_data,
                              interpolated_Y  = interpolated_Y,
                              X               = X,
                              list_indx_lst   = list_indx_lst,
                              filter_cs       = filter_cs,
                              outing_grid     = outing_grid,
                              cov_lev         = cov_lev,
                              filter.number   = filter.number,
                              family          = family,
                              ind_analysis    = ind_analysis,
                              post_processing = post_processing

    )
   multfsusie.obj$runtime <- proc.time()-pt
  return(multfsusie.obj)

}
