
#' @title Find affected region for effect l omics k
#' @param multfsusie.obj a multfsusie object
#' @param l the effect number
#' @param k a digit
#' @param \ldots other arguments
#' @export
#' @keywords internal

affected_reg_effect <- function(multfsusie.obj,l,k , ...)
  UseMethod("affected_reg_effect")


#' @rdname affected_reg_effect
#
#' @method affected_reg_effect multfsusie
#
#' @export affected_reg_effect.multfsusie
#
#' @export
#' @keywords internal
#' @importFrom stats complete.cases


affected_reg_effect.multfsusie <- function( multfsusie.obj, l,k , ... ){
  if(is.null(multfsusie.obj$fitted_wc)){
    return(multfsusie.obj)
  }

  outing_grid <- multfsusie.obj$outing_grid

  reg <-  list()
  h <- 1
  for (   l in 1:length(multfsusie.obj$cs)){

    pos_up <-  which(multfsusie.obj$cred_band[[l]][[k]][1,]<0)
    pos_low <- which(multfsusie.obj$cred_band[[l]][[k]][2,]>0)


    reg_up <- split( pos_up,cumsum(c(1,diff( pos_up)!=1)))

    reg_low <- split( pos_low,cumsum(c(1,diff( pos_low)!=1)))
    for( o in 1:length(reg_up)){
      reg[[h]] <- c(l, outing_grid[reg_up[[o]][1]], outing_grid[reg_up[[o]][length(reg_up[[o]])]])


      h <- h+1
    }
    for( k in 1:length(reg_low )){
      reg[[h]] <- c(l, outing_grid[reg_low [[o]][1]], outing_grid[reg_low [[o]][length(reg_low [[o]])]])


      h <- h+1
    }


  }


  reg <-  do.call(rbind, reg)
  colnames(reg) <- c("CS", "Start","End")
  reg <- as.data.frame(reg)
  reg <- reg[order(reg$CS, reg$Start),]
  reg <- reg[stats::complete.cases(reg),]
  return(reg)
}

#' @title Compute partial residuals
#' @param multfsusie.obj a multfsusie object
#' @param \ldots other arguments
#' @export
#' @keywords internal

cal_partial_resid <- function(multfsusie.obj,...)
  UseMethod("cal_partial_resid")



#' @rdname cal_partial_resid
#
#' @method cal_partial_resid multfsusie
#
#' @export cal_partial_resid.multfsusie
#
#' @export
#' @keywords internal


cal_partial_resid.multfsusie <- function(multfsusie.obj = multfsusie.obj,
                                         l = l,
                                         X = X,
                                         Y = Y,
                                         list_indx_lst  = list_indx_lst,...){

  L <- multfsusie.obj$L

  update_Y <- list(Y_u=NULL,
                   Y_f=NULL)

  if( !is.null(Y$Y_u)){
    if(L>1){
      id_L <- (1:L)[ - ( (l%%L)+1) ]#Computing residuals R_{l+1} by removing all the effect except effect l+1
      update_Y$Y_u   <- Y$Y_u - Reduce("+", lapply(id_L, function(k)
                                                    pred_partial_u(multfsusie.obj,k,X)
                                                    )
                                        )

    }else{
      id_L <- 1
      update_Y$Y_u   <- Y$Y_u - pred_partial_u(multfsusie.obj,1,X)
    }
  }
  if(!is.null(Y$Y_f)){
    update_Y$Y_f   <-   lapply(1:length(Y$Y_f),
                               function(k)
                                 cal_partial_resid_sub(multfsusie.obj,
                                                      l=l,
                                                      X=X,
                                                      D= Y$Y_f[[k]][,-list_indx_lst[[k]][[length(list_indx_lst[[k]])]]],
                                                      C= Y$Y_f[[k]][,list_indx_lst[[k]][[length(list_indx_lst[[k]])]]],
                                                      indx_lst = list_indx_lst[[k]],
                                                      cord=k
                                                      )
                               )

  }


  return(update_Y)

}


cal_partial_resid_sub <- function( multfsusie.obj, l, X, D, C, indx_lst,cord){
  L <- multfsusie.obj$L

  if(L>1){
    id_L <- (1:L)[ - ( (l%%L)+1) ]#Computing residuals R_{l+1} by removing all the effect except effect l+1


    update_D  <-  D - Reduce("+", lapply  ( id_L, function(l) (X*rep(multfsusie.obj$alpha[[l]],
                                                                     rep.int(dim(X)[1],dim(X)[2]))) %*% (multfsusie.obj$fitted_wc[[l]][[cord]][,-indx_lst[[length(indx_lst)]]])
                                            )
                             )
    update_C  <-  C - Reduce("+", lapply  ( id_L, function(l) (X*rep(multfsusie.obj$alpha[[l]],
                                                                     rep.int(dim(X)[1],dim(X)[2]))) %*% multfsusie.obj$fitted_wc[[l]][[cord]][,indx_lst[[length(indx_lst)]]]
                                            )
                             )
    update_Y  <- cbind(  update_D, update_C)
  }else{
    id_L <- 1


    update_D  <-  D - Reduce("+", lapply  ( id_L, function(l) (X*rep(multfsusie.obj$alpha[[l]], rep.int(dim(X)[1],dim(X)[2]))) %*% (multfsusie.obj$fitted_wc[[l]][[cord]][,-indx_lst[[length(indx_lst)]]])   ) )
    update_C  <-  C - Reduce("+", lapply  ( id_L, function(l) (X*rep(multfsusie.obj$alpha[[l]], rep.int(dim(X)[1],dim(X)[2]))) %*% multfsusie.obj$fitted_wc[[l]][[cord]][,indx_lst[[length(indx_lst)]]] ) )
    update_Y  <- cbind(  update_D, update_C)
  }
  return(update_Y)

}



#
#' @title Check purity credible sets
#
#' @param multfsusie.obj a multfsusie object defined by init_multfsusie_obj function
#' @param min_purity minimal purity within a CS
#' @param X matrix of covariate
#' @return a multfsusie.obj without "dummy" credible s
#
#' @export
#
#
check_cs <- function(multfsusie.obj, min_purity=0.5,X,...)
  UseMethod("check_cs")



#' @rdname check_cs
#
#' @method check_cs multfsusie
#
#' @export check_cs.multfsusie
#
#' @export
#' @keywords internal


check_cs.multfsusie <- function(multfsusie.obj, min_purity=0.5,X, ... )
{


  dummy.cs <- which_dummy_cs(multfsusie.obj, min_purity=min_purity,X)


  if( length(dummy.cs)==0)
  {
    return(multfsusie.obj)
  }else{
    if(length(dummy.cs)==multfsusie.obj$L) #avoid returning empty results
    {
      dummy.cs <- dummy.cs[-length(dummy.cs)]
    }
    multfsusie.obj <- discard_cs( multfsusie.obj,cs=dummy.cs, out_prep= TRUE)
    return(multfsusie.obj)
  }


}






create_dummy_susiF <- function( multfsusie.obj   ){



  G_prior <- list()
  class(G_prior) <- "mixture_normal_per_scale"
  susiF.obj <- list(cs=multfsusie.obj$cs ,
                    G_prior = G_prior,
                    alpha= multfsusie.obj$alpha,
                    csd_X= multfsusie.obj$csd_X
  )


  class(susiF.obj)  <- "susiF"

  return( susiF.obj)
}



#' @title Discard credible sets
#
#' @param multfsusie.obj a multfsusie object defined by init_multfsusie_obj function
#
#' @param cs vector of integer containing the credible sets to discard
#' @param ... Additional arguments passed to other functions.
#
#' @return a multfsusie.obj without "dummy" credible sets
#
#' @export
#' @keywords internal




discard_cs <- function(multfsusie.obj, cs,...)
  UseMethod("discard_cs")



#' @rdname discard_cs
#
#' @method discard_cs multfsusie
#
#' @export discard_cs.multfsusie
#
#' @export
#

discard_cs.multfsusie <- function(multfsusie.obj, cs, out_prep=FALSE, ...)
{

  if( length(cs)==multfsusie.obj$L){
    cs <- cs[-1]
    if(length(cs)==0){
      return(multfsusie.obj)
    }
  }
  if( length(cs)>0){
  multfsusie.obj$alpha       <-  multfsusie.obj$alpha[ -cs]
  multfsusie.obj$lBF         <-  multfsusie.obj$lBF[ -cs]
  if(!is.null(multfsusie.obj$G_prior$G_prior_f)){
    multfsusie.obj$fitted_wc   <-  multfsusie.obj$fitted_wc[ -cs]
    multfsusie.obj$fitted_wc2  <-  multfsusie.obj$fitted_wc2[ -cs]
  }
  if(!is.null(multfsusie.obj$G_prior$G_prior_u)){
    multfsusie.obj$fitted_u  <-  multfsusie.obj$fitted_u[-cs]
    multfsusie.obj$fitted_u2 <-  multfsusie.obj$fitted_u2[-cs]
  }


  if(out_prep){
    #multfsusie.obj$fitted_func <-  multfsusie.obj$fitted_func[ -cs]
  }else{
    multfsusie.obj$greedy_backfit_update <- TRUE
    multfsusie.obj$KL                    <- multfsusie.obj$KL[ -cs]
    multfsusie.obj$ELBO                  <- -Inf
  }

  multfsusie.obj$cs          <-  multfsusie.obj$cs[ -cs]

  multfsusie.obj$est_pi           <-  multfsusie.obj$est_pi[ -cs]
  multfsusie.obj$L                <-  multfsusie.obj$L -length(cs)
}

  return(multfsusie.obj)
}




# @title Expand multfsusie.obj by adding L_extra effect
#
# @param multfsusie.obj a multfsusie.obj
#
# @param L_extra numeric a number of effect to add
#
# @return a multfsusie.obj a L_extra effect. Note the the number of effect of the multfsusie.obj cannot exceed the number the user upper bound
expand_multfsusie_obj <- function(multfsusie.obj,L_extra)
{

  L_extra <- ifelse (  multfsusie.obj$L_max - (multfsusie.obj$L+L_extra) <0 ,#check if we are adding more effect that maximum specified by user
                       abs(multfsusie.obj$L_max -(multfsusie.obj$L)),
                       L_extra

                    )

  if( L_extra==0){
    return(multfsusie.obj)
  }else{
    L_old <- multfsusie.obj$L
    L_new <- multfsusie.obj$L+L_extra
    multfsusie.obj$L <- ifelse(L_new<(multfsusie.obj$P+1),L_new,multfsusie.obj$P)
    l=1
    k=1

    for ( l in (L_old+1):multfsusie.obj$L )
    {

      if( !is.null(multfsusie.obj$fitted_wc)){

        multfsusie.obj$fitted_wc[[l]]  <-    lapply( 1:length(multfsusie.obj$n_wac), function(j)
                                                      matrix( 0,
                                                              ncol = multfsusie.obj$n_wac[[j]] ,
                                                              nrow =  multfsusie.obj$P)
                                                     )
        multfsusie.obj$fitted_wc2[[l]] <-    lapply( 1:length(multfsusie.obj$n_wac), function(j)
                                                       matrix( 1,
                                                              ncol = multfsusie.obj$n_wac[[j]],
                                                              nrow = multfsusie.obj$P)
                                                     )

         }
      if(!is.null(multfsusie.obj$fitted_u)){
        multfsusie.obj$fitted_u[[l]]        <- 0*multfsusie.obj$fitted_u[[1]]
        multfsusie.obj$fitted_u2[[l]]       <- 0*multfsusie.obj$fitted_u2[[1]]

      }

      multfsusie.obj$alpha [[l]]           <- rep(0, length(multfsusie.obj$alpha [[1]]))
      multfsusie.obj$cs[[l]]               <- list()
      multfsusie.obj$est_pi [[l]]          <- multfsusie.obj$est_pi[[1]]
      multfsusie.obj$lBF[[l]]              <- rep(NA, length( multfsusie.obj$lBF[[1]]))
      multfsusie.obj$KL                    <- rep(NA,multfsusie.obj$L)
      multfsusie.obj$ELBO                  <- c()
    }
    multfsusie.obj$n_expand <- multfsusie.obj$n_expand+1
    multfsusie.obj$greedy_backfit_update <- TRUE

    if( multfsusie.obj$L== multfsusie.obj$L_max){
      multfsusie.obj$greedy=FALSE
    }
    return(multfsusie.obj)
  }

}

#'@title Compute refined estimate using HMM regression
#'
#' @description    Compute refined estimate using HMM regression
#' @param multfsusie.obj  a multfsusie.obj object
#'
#' @param Y  multfsusie.obj data object
#'
#' @param X matrix containing the covariates
#'
#' @param verbose logical
#'
#' @param maxit maximum number of iterations
#'
#' @param ind_analysis output of which_notNA_pos
#'
#' @param ... Additional arguments passed to other functions.
#'
#' @export


HMM_regression<- function (multfsusie.obj,Y,X,verbose, maxit,...)
  UseMethod("HMM_regression")


#' @rdname HMM_regression
#'
#' @method HMM_regression multfsusie
#'
#' @export HMM_regression.multfsusie
#'
#' @export
#'


HMM_regression.multfsusie <- function(multfsusie.obj,Y,X ,verbose=TRUE,maxit=5, ind_analysis ,...   ){

  if(is.null(multfsusie.obj$fitted_wc)){
    return(multfsusie.obj)
  }
  if(verbose){
    print( "Fine mapping done, refining effect estimates using HMM regression")
  }

 #  browser()
  #renaming lfsr output
  multfsusie.obj$lfsr <- vector("list", length(multfsusie.obj$cs))

  for (l in seq_along(multfsusie.obj$lfsr)) {
    multfsusie.obj$lfsr[[l]] <- list(
      est_lfsr_functional = list()
    )
  }

  tl <- list()

  for ( k in 1: length(Y$Y_f)){
    tl[[k]] <- rep (0, ncol (Y$Y_f[[k]]))


  }
  tp <- list()
  for ( l in 1:length(multfsusie.obj$cs)){
    multfsusie.obj$fitted_func[[l]] <- tl

    tp[[l]] <- list()
    for ( k in 1: length(Y$Y_f)){
      tp[[l]][[k]]  <- rbind ( rep( 0 , ncol( Y$Y_f[[k]])),
                               rep( 0 , ncol( Y$Y_f[[k]])))
    }
  }
  multfsusie.obj$cred_band <- tp
  dummy_susiF.obj <- create_dummy_susiF(multfsusie.obj)

  for ( k in 1: length(Y$Y_f)){


    susiF.obj <- fsusieR::HMM_regression.susiF(    obj            = dummy_susiF.obj,
                                                    Y             = Y$Y_f[[k]][ind_analysis$idx_f[[k]],],
                                                    X             = X[ind_analysis$idx_f[[k]],, drop=FALSE],
                                                    verbose       = FALSE ,
                                                    fit_indval    = FALSE)

    for ( l in 1:length(multfsusie.obj$cs)){

      multfsusie.obj$fitted_func[[l]][[k]] <- susiF.obj$fitted_func[[l]]
      multfsusie.obj$lfsr[[l]]$est_lfsr_functional[[k]] <-  susiF.obj$lfsr_func[[l]]

    }


  }


  return(multfsusie.obj)
}


# @title Initialize a multsusie object
#
# @param L_max upper bound on the number of non zero coefficients An L-vector containing the indices of the
#   nonzero coefficients.
#
#
# @param G_prior prior object defined by init_prior_multsusie function
#
# @param Y list of matrices of outcomes
#
# @param X Matrix of covariates
#
# @param L_start number of effect to start with
# @pama type_mark an object generated by \code{\link{is.functional}} function
#
# @return A list with the following elements
# \item{fitted_wc}{ list of length L, each element contains the fitted wavelet coefficients of effect l}
# \item{fitted_wc2}{list of length L, each element contains the variance of the fitted wavelet coefficients of effect l}
# \item{alpha_hist}{ history of the fitted alpha value}
# \item{N}{ number of indidivual in the study}
# \item{sigma2}{residual variance}
# \item{n_wac}{number of wavelet coefficients}
# \item{ind_fitted_func}{fitted curves of each individual }
# \item{cs}{credible set}
# \item{pip}{Posterior inclusion probabilites}
# \item{G_prior}{a G_prior of the same class as the input G_prior, used for internal calculation}
# \item{lBF}{ log Bayes factor for the different effect}
# \item{KL}{ the KL divergence for the different effect}
# \item{ELBO}{ The evidence lower bound}
# \item{lfsr_wc}{Local fasle sign rate of the fitted wavelet coefficients}
# @export
init_multfsusie_obj <- function(L_max,
                                G_prior,
                                Y,
                                X,
                                type_mark,
                                L_start,
                                greedy,
                                backfit,
                                ind_analysis,
                                tol_null_prior,
                                lbf_min=0.1)
{
  sigma2          <- list()
  if(!is.null(Y$Y_f)){
  fitted_wc       <- list()
  fitted_wc2      <- list()
  n_wac           <- lapply(lapply(Y$Y_f,dim) ,`[[`, 2)
  sigma2$sd_f     <- rep( 1, length(n_wac ))
  lfsr_wc         <- list()
  }else {
    fitted_wc     <- NULL
    fitted_wc2    <- NULL
    n_wac         <- NULL
    sigma2$sd_f   <- NULL
    lfsr_wc       <- NULL
  }
  if(!is.null(Y$Y_u)){
    fitted_u    <- list()
    fitted_u2   <- list()
    lfsr_u      <- list()
    n_univ      <- ncol(Y$Y_u)


  }else{
    fitted_u    <- NULL
    fitted_u2   <- NULL
    sigma2$sd_u <- NULL
    lfsr_u      <- NULL
    n_univ      <- NULL
  }
  alpha           <- list()
  alpha_hist      <- list()
  ind_fitted_val  <- list()
  cs              <- list()
  pip             <- rep(0, dim(X)[2])
  est_pi          <- list()
  est_sd          <- list()
  L_max           <- L_max
  L               <- L_start
  G_prior         <- G_prior
  N               <- nrow(X)[1]
  n_cond          <- type_mark$ncond
  P               <- ncol(X)
  lBF             <- list()
  lBF_per_trait   <- list()
  KL              <- rep(NA,L)
  ELBO            <- c()
  P               <- ncol(X)
  mean_X          <- attr(X, "scaled:center")
  csd_X           <- attr(X, "scaled:scale")
  n_expand        <- 0 #number of greedy expansion
  greedy          <- greedy
  backfit         <- backfit
  greedy_backfit_update <- FALSE
  ind_analysis    <- ind_analysis
  for ( l in 1:L )
  {

    if(!is.null(Y$Y_f)){
      fitted_wc [[l]]       <-    lapply( 1:length(n_wac),
                                          function(j) matrix( 0,ncol= n_wac[[j]],
                                                              nrow = ncol(X))
                                          )
      fitted_wc2[[l]]       <-    lapply( 1:length(n_wac),
                                          function(j) matrix(  1,ncol= n_wac[[j]],
                                                               nrow = ncol(X))
                                          )
      lfsr_wc   [[l]]       <-    lapply(lapply(Y$Y_f,ncol) ,
                                         function(i) rep(1,i)
                                         )

    }
    if(!is.null(Y$Y_u)){
      fitted_u [[l]]       <-     matrix(0, ncol= ncol(Y$Y_u), nrow = ncol(X))
      fitted_u2[[l]]       <-     matrix(0, ncol= ncol(Y$Y_u), nrow = ncol(X))
      lfsr_u   [[l]]       <-     rep(1,ncol(Y$Y_u))
      sigma2$sd_u          <-     rep( 1, ncol(Y$Y_u))
    }


    alpha [[l]]          <-  rep(0, dim(X)[2])
    cs    [[l]]          <-  list()
    est_pi[[l]]          <-  get_pi_G_prior(G_prior)
    lBF   [[l]]          <-  rep(NA, ncol(X))

  }

 sigma2 <-  init_var_multf(Y)

  obj <- list( fitted_wc       = fitted_wc,
               fitted_wc2      = fitted_wc2,
               fitted_u        = fitted_u,
               fitted_u2       = fitted_u2,
               lBF             = lBF,
               lBF_per_trait   = lBF_per_trait,
               KL              = KL,
               ELBO            = ELBO,
               ind_fitted_val  = ind_fitted_val,
               G_prior         = G_prior,
               alpha_hist      = alpha_hist,
               N               = N,
               n_wac           = n_wac,
               n_univ          = n_univ,
               sigma2          = sigma2,
               P               = P,
               alpha           = alpha,
               cs              = cs,
               pip             = pip,
               est_pi          = est_pi,
               est_sd          = est_sd,
               mean_X          = mean_X,
               csd_X           = csd_X,
               L               = L,
               P               = P,
               L_max           = L_max,
               n_expand        = n_expand,
               greedy          = greedy,
               backfit         = backfit,
               greedy_backfit_update=greedy_backfit_update,
               ind_analysis    = ind_analysis,
               tol_null_prior  = tol_null_prior,
               lbf_min         =  lbf_min )

  class(obj) <- "multfsusie"
  return(obj)
}




#' @title Acces mixture proportion  for multfsusie.obj
#' @description see title
#' @param multfsusie.obj a multfsusie object
#' @param l effect to be accesed
#' @export
#' @keywords internal
get_pi <- function(multfsusie.obj,l, ...)
  UseMethod("get_pi")



#' @rdname get_pi
#
#' @method get_pi multfsusie
#
#' @export get_pi.multfsusie
#
#' @export
#' @keywords internal
get_pi.multfsusie <- function(multfsusie.obj, l, ...)
{

  if( l >  length(multfsusie.obj$est_pi))
  {
    stop("Error trying to access mixture proportion")
  }
  if( l < 1)
  {
    stop("Error l should be larger ")
  }
  out <- multfsusie.obj$est_pi[[l]]
  return(out)
}





#' @title Acces mixture proportion for null component for multfsusie.obj
#' @description see title
#' @param multfsusie.obj a multfsusie object
#' @param l effect to be accesed
#' @export
#' @keywords internal

get_pi0<- function(multfsusie.obj,l, ...)
  UseMethod("get_pi0")


#' @rdname get_pi0
#
#' @method get_pi0 multfsusie
#
#' @export get_pi0.multfsusie
#' @keywords internal




get_pi0.multfsusie <-function(multfsusie.obj, l, ... ){

  if (missing( l)){


    if( is.null(multfsusie.obj$G_prior$G_prior_u))
    {
      out <- list()
      for ( l in 1:1:multfsusie.obj$L ){

        if (is.null( multfsusie.obj$n_univ)){
          pi0_u <-NULL
        }else{
          pi0_u <- lapply( 1: multfsusie.obj$n_univ ,
                           function(k)
                             multfsusie.obj$est_pi[[l]]$est_pi_u[[k]]
          )
        }

        if (is.null( multfsusie.obj$n_wac)){
          pi0_f <-NULL
        }else{
          pi0_f <- lapply( 1:length(multfsusie.obj$n_wac),
                           function(k)
                             sapply(multfsusie.obj$est_pi[[l]]$est_pi_f[[k]],"[[",1)
          )
        }

        out[[l]]  <- list( pi0_u=pi0_u,
                           pi0_f = pi0_f)
    }



    }
  }else{

    if (is.null( multfsusie.obj$n_univ)){
      pi0_u <-NULL
    }else{
      pi0_u <- lapply( 1: multfsusie.obj$n_univ ,
                       function(k)
                         multfsusie.obj$est_pi[[l]]$est_pi_u[[k]]
      )
    }

    if (is.null( multfsusie.obj$n_wac)){
      pi0_f <-NULL
    }else{
      pi0_f <- lapply( 1:length(multfsusie.obj$n_wac),
                       function(k)
                         sapply(multfsusie.obj$est_pi[[l]]$est_pi_f[[k]],"[[",1)
      )
    }



      out   <- list( pi0_u=pi0_u,
                         pi0_f = pi0_f)
    }

  return(out)
}




#' @title Acces prior object for multfsusie.obj
#' @description see title
#' @param multfsusie.obj a multfsusie object
#'  @param ... Additional arguments passed to other functions.
#' @export
#' @keywords internal


get_G_prior<- function(multfsusie.obj, ...)
  UseMethod("get_G_prior")


#' @rdname get_G_prior
#
#' @method get_G_prior multfsusie
#
#' @export get_G_prior.multfsusie
#
#' @export
#' @keywords internal
get_G_prior.multfsusie <- function(multfsusie.obj, ...){
  out <- multfsusie.obj$G_prior
  return(out)
}




#' @title Acces log Bayes factors for multfsusie.obj
#' @description see title
#' @param multfsusie.obj a multfsusie object
#' @param l effect to be accesed
#' @export
#' @keywords internal


get_lBF<- function(multfsusie.obj,l, ...)
  UseMethod("get_lBF")
#' @rdname get_lBF
#
#' @method get_lBF multfsusie
#
#' @export get_lBF.multfsusie
#
#' @export
#' @keywords internal
get_lBF.multfsusie <- function(multfsusie.obj,l, ...){
  out <- multfsusie.obj$lBF[[l]]
  return(out)
}



# @title Compute posterior mean of the fitted effect
#
# @param G_prior a multfsusie_prior object
#
# @param effect_estimate an object generated by \link{\code{cal_Bhat_Shat_multfsusie }}
# @param list_indx_lst list of list generated by \code{\link{gen_wavelet_indx}} for the given level of resolution
# @param lowc_wc
#
#
#
# @return  an object of the same form as effect_estimate, which corresponds to the posterior mean.
get_post_effect_multfsusie <- function(G_prior, effect_estimate,lBF, list_indx_lst=NULL, low_trait = NULL, e=0.001){

  out <- list( res_u = NULL,
               res_f  = NULL)
 if( !is.null(effect_estimate$res_u)){


     out$res_u$Bhat  <- do.call(
                                  cbind,
                                  lapply( 1:length(G_prior$G_prior_u),
                                          function( k) get_post_mean_u(G_prior$G_prior_u[[k]],
                                                                       effect_estimate$res_u$Bhat[,k],
                                                                       effect_estimate$res_u$Shat[,k],
                                                                       low_u = ifelse(k%in%low_trait$low_u, TRUE,FALSE)
                                                                        )
                                          )
                                  )
     out$res_u$Shat  <- do.call(
                                  cbind,
                                  lapply( 1:length(G_prior$G_prior_u),
                                          function( k) get_post_sd_u(G_prior$G_prior_u[[k]],
                                                                     effect_estimate$res_u$Bhat[,k],
                                                                     effect_estimate$res_u$Shat[,k],
                                                                     low_u = ifelse(k%in%low_trait$low_u, TRUE,FALSE)
                                                                        )
                                         )
                                  )

   }

  if( !is.null(effect_estimate$res_f)){
   out$res_f <- lapply( 1:length(G_prior$G_prior_f), function(k)list_post_mean_sd(G_prior$G_prior_f[[k]],
                                                               effect_estimate$res_f[[k]]$Bhat,
                                                               effect_estimate$res_f[[k]]$Shat,
                                                               indx_lst =list_indx_lst[[k]],
                                                               lBF=lBF,
                                                               lowc_wc= low_trait$low_wc[[k]],
                                                               e = e)
                        )
  }
  return( out)
}





get_post_F  <- function(multfsusie.obj,l, ...)
  UseMethod("get_post_F")



# @rdname get_post_F
#
# @method get_post_F multfsusie
#
# @export get_post_F.multfsusie
#
#' @export
#' @keywords internal
get_post_F.multfsusie <- function (multfsusie.obj,l, ... ){

  if(missing(l))
  {
    if(!is.null(multfsusie.obj$fitted_wc)){
      out_f <-  lapply(1: length(multfsusie.obj$n_wac) ,
                       function(k) Reduce("+",
                                          lapply(1:multfsusie.obj$L,
                                                 FUN=function(l) multfsusie.obj$alpha[[l]] * multfsusie.obj$fitted_wc[[l]][[k]]
                                          )
                       )
      )
    }else{
      out_f <- NULL
    }
    if(!is.null(multfsusie.obj$fitted_u)){
      out_u <-  Reduce("+",
                       lapply(1:multfsusie.obj$L,
                              function(l)
                                multfsusie.obj$alpha[[l]] * multfsusie.obj$fitted_u[[l]]
                       )
      )
    }else{
      out_u <- NULL
    }

  }else{

      if(!is.null(multfsusie.obj$fitted_wc)){
              out_f <-  lapply(1: length(multfsusie.obj$n_wac),
                               function (k) multfsusie.obj$alpha[[l]] * multfsusie.obj$fitted_wc[[l]][[k]]
                               )
      }else{
              out_f <- NULL
      }
      if(!is.null(multfsusie.obj$fitted_u)){
             out_u <-  multfsusie.obj$alpha[[l]] * multfsusie.obj$fitted_u[[l]]
      }else{
             out_u <- NULL
      }
  }

  out_post <- list(post_u = out_u,
                   post_f   = out_f
                   )
  return(out_post)
}







get_post_F2  <- function(multfsusie.obj,l, ...)
  UseMethod("get_post_F2")






# @rdname get_post_F2
#
# @method get_post_F2 multfsusie
#
# @export get_post_F2.multfsusie
#
# @export
#' @export
#' @keywords internal
get_post_F2.multfsusie <- function (multfsusie.obj,l, ...){

  if(missing(l))
  {
    if(!is.null(multfsusie.obj$fitted_wc)){
      out_f <-  lapply(1: length(multfsusie.obj$n_wac) ,
                       function(k) Reduce("+",
                                          lapply(1:multfsusie.obj$L,
                                                 FUN=function(l) multfsusie.obj$alpha[[l]] *(multfsusie.obj$fitted_wc[[l]][[k]]^2+ multfsusie.obj$fitted_wc2[[l]][[k]])
                                          )
                       )
      )
    }else{
      out_f <- NULL
    }
    if(!is.null(multfsusie.obj$fitted_u)){
      out_u <-  Reduce("+",
                       lapply(1:multfsusie.obj$L,
                              function(l)
                                multfsusie.obj$alpha[[l]] * (multfsusie.obj$fitted_u2[[l]]+multfsusie.obj$fitted_u[[l]]^2    )
                       )
      )
    }else{
      out_u <- NULL
    }

  }else{

    if(!is.null(multfsusie.obj$fitted_wc)){
      out_f <-  lapply(1: length(multfsusie.obj$n_wac),
                       function (k)  multfsusie.obj$alpha[[l]] *(multfsusie.obj$fitted_wc[[l]][[k]]^2+ multfsusie.obj$fitted_wc2[[l]][[k]])
                       )


    }else{
      out_f <- NULL
    }
    if(!is.null(multfsusie.obj$fitted_u)){
      out_u <-   multfsusie.obj$alpha[[l]] * (multfsusie.obj$fitted_u2[[l]]+multfsusie.obj$fitted_u[[l]]^2    )

    }else{
      out_u <- NULL
    }
  }

  out_post <- list(post_u_sd2 = out_u,
                   post_f_sd2    = out_f
  )
  return(out_post)
}


#' @title Compute expected sum of squares
#' @description  See title
#' @param multfsusie.obj a multfsusie object
#' @param Y matrix of observation
#' @param X matri of covaraites
#' @export
#' @keywords internal
get_ER2  <- function(multfsusie.obj,Y,X, ...)
  UseMethod("get_ER2")


#' @rdname get_ER2
#
#' @method get_ER2 multfsusie
#
#' @export get_ER2.multfsusie
#
#' @export
#' @keywords internal

get_ER2.multfsusie = function (  multfsusie.obj,Y, X,ind_analysis, ... ) {
  postF <- get_post_F(multfsusie.obj )# J by N matrix
  #Xr_L = t(X%*% postF)
  postF2 <- get_post_F2(multfsusie.obj ) # Posterior second moment.

  ER2 <-  list()
  if(! is.null(Y$Y_u)){

    if( missing(ind_analysis)){
      ER2$uni <-  do.call( c,
                           lapply(1:ncol( Y$Y_u),
                                  function(k) sum((Y$Y_u[,k] - X%*%postF$post_u[,k] )^2) -sum( postF$post_u_sd2[,k]^2) +sum( postF2$post_u_sd2[,k])
                           )
      )
    }else{
      ER2$uni <-  do.call( c,
                           lapply(1:ncol( Y$Y_u),
                                  function(k) sum((Y$Y_u[ind_analysis$idx_u[[k]],k] - X[ind_analysis$idx_u[[k]],,drop=FALSE]%*%postF$post_u[,k] )^2) -sum( postF$post_u_sd2[,k]^2) +sum( postF2$post_u_sd2[,k])
                           )
      )
    }




  }else
  {
    ER2$uni <- NULL
  }
  if( !is.null(Y$Y_f))
  {
    if( missing(ind_analysis)){
      ER2$f <-  do.call( c,
                         lapply(1:length( Y$Y_f),
                                function(k) sum((Y$Y_f[[k]] - X%*%postF$post_f[[k]])^2)  -sum(postF$post_f[[k]]^2) + sum( postF2$post_f_sd2 [[k]])
                         )
      )
    }else{
      ER2$f <-  do.call( c,
                         lapply(1:length( Y$Y_f),
                                function(k) sum((Y$Y_f[[k]][ind_analysis$idx_f[[k]],] - X[ind_analysis$idx_f[[k]],, drop=FALSE]%*%postF$post_f[[k]])^2)  -sum(postF$post_f[[k]]^2) + sum( postF2$post_f_sd2 [[k]])
                         )
      )
    }


  }else{
    ER2$f <- NULL
  }
  return(ER2)
}




#' @title Update  susiF via greedy search or backfit
#
#' @param multfsusie.obj a susiF object defined by init_multfsusie_obj function
#
#' @return multfsusie.obj object
#
#' @importFrom fsusieR cal_cor_cs
#' @export get_pi0.multfsusie
#' @keywords internal
#
greedy_backfit  <-  function(multfsusie.obj, verbose,cov_lev,X,min_purity, ...  )
  UseMethod("greedy_backfit")

#' @rdname greedy_backfit
#'
#' @method greedy_backfit multfsusie
#'
#' @export greedy_backfit.multfsusie
#
#' @export get_pi0.multfsusie
#' @keywords internal

greedy_backfit.multfsusie <-  function(multfsusie.obj,verbose,cov_lev,X,min_purity, ...  )
{


  multfsusie.obj <- update_alpha_hist(multfsusie.obj)
  if(!(multfsusie.obj$greedy)&!(multfsusie.obj$backfit))
  {
    return(multfsusie.obj)
  }
  multfsusie.obj <- update_cal_cs(multfsusie.obj,
                                  cov_lev=cov_lev)

  dummy.cs <-  which_dummy_cs(multfsusie.obj,
                              min_purity = min_purity,
                              median_crit=TRUE,
                              X=X,
                              lbf_min=multfsusie.obj$lbf_min)
  if(multfsusie.obj$backfit & (length(dummy.cs)>0)){

    multfsusie.obj$greedy <- FALSE
    if(length(dummy.cs)== multfsusie.obj$L){###  avoid returning empty object
      dummy.cs <- dummy.cs[-1]
      multfsusie.obj$backfit <- FALSE
    }
    if( length(dummy.cs)==0  )
    {
      multfsusie.obj$backfit <- FALSE
    }else{
      temp_L <- multfsusie.obj$L


      multfsusie.obj <- discard_cs(multfsusie.obj,
                                   cs= dummy.cs,
                                   out_prep= FALSE
      )

      merge_effect (multfsusie.obj , verbose = verbose)
      if(verbose){
        print( paste( "Discarding ",(temp_L- multfsusie.obj$L), " effects"))
      }
    }
    return(multfsusie.obj)

  }##Conditions for stopping greedy search
  if(  (multfsusie.obj$L>multfsusie.obj$L_max))
  {

    multfsusie.obj$greedy <- FALSE

    merge_effect (multfsusie.obj , verbose = verbose)

    multfsusie.obj <- discard_cs(multfsusie.obj,
                                 cs= (multfsusie.obj$L_max+1):multfsusie.obj$L,
                                 out_prep= FALSE
    )




    multfsusie.obj  <- merge_effect (multfsusie.obj , verbose = verbose)
    if(verbose){
      print( paste( "Discarding ",(multfsusie.obj$L_max- multfsusie.obj$L), " effects"))
      print( "Greedy search and backfitting done")
    }

  }

  if( length(dummy.cs)==0& !( multfsusie.obj$greedy))#### TO double
  {

    multfsusie.obj$backfit <- FALSE
  }

  if(!(multfsusie.obj$greedy )&!(multfsusie.obj$backfit ) ){

    multfsusie.obj  <- merge_effect (multfsusie.obj , verbose = verbose)

    if(verbose){
      print( paste( "Discarding ",(multfsusie.obj$L_max- multfsusie.obj$L), " effects"))
      print( "Greedy search and backfitting done")
    }
    multfsusie.obj <- update_alpha_hist(multfsusie.obj,discard = TRUE)
    multfsusie.obj$greedy_backfit_update <- FALSE

    return(multfsusie.obj)
  }
  if(multfsusie.obj$greedy & (length(dummy.cs)==0)){

    tt <- multfsusie.obj$L_max -multfsusie.obj$L

    temp <- min( ifelse(tt>0,tt,0 ) , 7)

    if(temp==0){
      multfsusie.obj  <- merge_effect (multfsusie.obj , verbose = verbose)

      if(verbose){

        print( paste( "Discarding ",(multfsusie.obj$L_max- multfsusie.obj$L), " effects"))
        print( "Greedy search and backfitting done")
      }
      multfsusie.obj <- update_alpha_hist(multfsusie.obj,discard = TRUE)
      multfsusie.obj$greedy_backfit_update <- FALSE
      multfsusie.obj$backfit <- FALSE
      multfsusie.obj$greedy <- FALSE
      return(multfsusie.obj)
    }


    if(verbose){
      print( paste( "Adding ", temp, " extra effects"))
    }



    multfsusie.obj  <- merge_effect (multfsusie.obj , verbose = verbose)



    multfsusie.obj <- expand_multfsusie_obj(multfsusie.obj,L_extra = 7)
    return(multfsusie.obj)
  }

}



# @title formatting function for output of posterior quantities
# @description formatting function for output of posterior quantities
# @param G_prior a mixutre_per_scale prior
# @param  Bhat matrix of estimated mean
# @param Shat matrix of estimated sd
# @importFrom fsusieR post_mat_mean
# @importFrom fsusieR post_mat_sd
#

list_post_mean_sd <- function(G_prior, Bhat,Shat,lBF,  indx_lst, lowc_wc=NULL,e=0.001)
{
  out <- list (Bhat= fsusieR::post_mat_mean( G_prior ,
                                                 Bhat,
                                                 Shat,
                                                 lBF      = lBF,
                                                 indx_lst = indx_lst,
                                                 lowc_wc  = lowc_wc,
                                                 e=e),
               Shat=fsusieR:: post_mat_sd(   G_prior ,
                                                 Bhat,
                                                 Shat,
                                                 lBF      = lBF,
                                                 indx_lst = indx_lst,
                                                 lowc_wc  = lowc_wc,
                                                 e=e)
  )
  return(out)

}



#' @title Merging effect function
#
#' @param multfsusie.obj a susiF object defined by init_multfsusie_obj function
#
#' @param  verbose
#
#
#
#' @return  a multfsusie  object
#' @export
#' @keywords internal
merge_effect <- function( multfsusie.obj,  verbose , ...)
  UseMethod("merge_effect")

#' @rdname merge_effect
#
#' @method merge_effect multfsusie
#
#' @export merge_effect.multfsusie
#
#' @export
#' @keywords internal

merge_effect.multfsusie  <-  function(multfsusie.obj, verbose = FALSE) {

  if (multfsusie.obj$L < 2) return(multfsusie.obj)

  to_drop <- integer(0)

  for (i in 1:(multfsusie.obj$L - 1)) {
    for (j in (i + 1):multfsusie.obj$L) {

      if (!(i %in% to_drop || j %in% to_drop)){
        rel <-  fsusieR::cs_relation(multfsusie.obj$cs[[i]],
                                     multfsusie.obj$cs[[j]])

        if(! (rel == "none")){
          largest <- fsusieR::which_cs_largeBF(multfsusie.obj, i, j)
          smallest  <- setdiff(c(i, j), largest)

          to_drop <- c(to_drop,  smallest )

          if (verbose) {
            message(
              sprintf(
                "Trivial merge: dropping effect %d (CS %s effect %d)",
                smallest,
                ifelse(rel == "identical", "identical to", "nested in"),
                largest
              )
            )
          }
        }


      }


    }
  }

  to_drop <- sort(unique(to_drop))
  if (length(to_drop) == 0) return(multfsusie.obj)

  discard_cs(multfsusie.obj, cs = to_drop, out_prep = FALSE)
}







# @title Updates CS names for output
#
# @param multfsusie.obj a multfsusie object defined by init_multfsusie_obj function
#
# @param X matrix of size N by p

name_cs <- function(multfsusie.obj,X,...)
  UseMethod("name_cs")

# @rdname name_cs
#
# @method name_cs multfsusie
#
# @export name_cs.multfsusie
#
#' @export
#' @keywords internal

name_cs.multfsusie <- function(multfsusie.obj,X,...){

  if( length(colnames(X))==ncol(X)){

    for (l in 1: length(multfsusie.obj$cs)){
      names(multfsusie.obj$cs[[l]]) <- colnames(X)[multfsusie.obj$cs[[l]]]
    }

  }
  return(multfsusie.obj)
}




#' @title Preparing output of main multfsusie function
#
#' @param multfsusie.obj a multfsusie object defined by init_multfsusie_obj function
#
#' @param Y  data
#' @param interpolated_Y interpolated functional data
#' @param X matrix of size N by p
#
#' @param list_indx_lst list generated by gen_wavelet_indx for the given level of resolution
#
#' @param filter_cs logical, if TRUE filter the credible set (removing low purity cs and cs with estimated prior equal to 0)
#'@param post_processing the chosen postprocessing
#'@param verbose  true or false
#' @export
#' @keywords internal
out_prep <- function(multfsusie.obj,
                     Y,
                     interpolated_Y,
                     X,
                     list_indx_lst,
                     filter_cs,
                     outing_grid ,
                     cov_lev ,
                     filter.number ,
                     family ,
                     ind_analysis ,
                     post_processing,...)
UseMethod("out_prep")

#' @rdname out_prep
#
#' @method out_prep multfsusie
#
#' @export out_prep.multfsusie
#
#' @export
#' @keywords internal



out_prep.multfsusie <- function(multfsusie.obj,
                                Y,
                                interpolated_Y,
                                X,
                                list_indx_lst,
                                filter_cs,
                                outing_grid,
                                cov_lev,
                                filter.number = 10,
                                family = "DaubLeAsymm",
                                ind_analysis ,
                                post_processing="smash",
                                verbose=TRUE,
                                ... )
{
  multfsusie.obj <-  update_cal_pip(multfsusie.obj)
  multfsusie.obj <-  update_cal_fit_u(multfsusie.obj )

  if(post_processing== "none"){
    multfsusie.obj <-  update_cal_fit_func(multfsusie.obj,list_indx_lst)

  }
  if( post_processing== "smash"){

    multfsusie.obj <-  smash_regression(multfsusie.obj = multfsusie.obj,
                                     Y              = interpolated_Y,
                                     X              = X,
                                     verbose        = verbose ,
                                     ind_analysis   = ind_analysis,
                                     filter.number  = filter.number,
                                     family         = family)
  }
  if( post_processing== "TI"){
    print("TIIII")
    multfsusie.obj <-  TI_regression(multfsusie.obj = multfsusie.obj,
                                     Y              = interpolated_Y,
                                     X              = X,
                                     verbose        = verbose ,
                                     ind_analysis   = ind_analysis,
                                     filter.number  = filter.number,
                                     family         = family)
  }
  if(post_processing=="HMM") {
    multfsusie.obj <-  HMM_regression(multfsusie.obj = multfsusie.obj,
                                      Y               = interpolated_Y,
                                      ind_analysis    = ind_analysis,
                                      X               = X,
                                      verbose         = verbose )
  }





  multfsusie.obj <- update_cal_cs(multfsusie.obj,
                                  cov_lev=cov_lev)

  multfsusie.obj <-  name_cs(multfsusie.obj,X)
  if(filter_cs)
  {
    multfsusie.obj <- check_cs(multfsusie.obj,min_purity=0.5,X=X)
  }


  if(!is.null( Y$Y_u)){

    multfsusie.obj$outing_grid <- outing_grid
  }
  multfsusie.obj$purity      <- fsusieR::cal_purity(l_cs= multfsusie.obj$cs, X=X)


  return( multfsusie.obj)
}


pred_partial_u <- function( multfsusie.obj, l, X )
{
  tbeta <-multfsusie.obj$alpha[[l]] *multfsusie.obj$fitted_u[[l]]
  pred_l   <- X%*% ( tbeta)
  return(pred_l)
}






#' @title Update alpha   susiF mixture proportion of effect l
#
#' @param multfsusie.obj a multfsusie  object defined by init_multfsusie_obj function
#
#' @param l integer larger or equal to 1. Corresponds to the effect to be accessed
#
#' @param alpha  vector of p alpha values summing up to one
#
#' @return multfsusie.obj object
#
#' @export
#' @keywords internal

update_alpha  <-  function(multfsusie.obj, l, alpha,... )
  UseMethod("update_alpha")

#' @rdname update_alpha
#
#' @method update_alpha multfsusie
#
#' @export update_alpha.multfsusie
#
#' @export
#' @keywords internal
update_alpha.multfsusie <-  function(multfsusie.obj, l, alpha, ... )
{
  multfsusie.obj$alpha[[l]] <- alpha

  return( multfsusie.obj)
}



#' @title Update alpha_hist   multfsusie object
#
#' @param multfsusie.obj a multfsusie object defined by init_multfsusie_obj function
#
#' @param  discard logical set to FALSE by default, if true remove element of history longer than L
#
#' @return multfsusie object
#
#' @export
#' @keywords internal

update_alpha_hist  <-  function(multfsusie.obj, discard, ... )
  UseMethod("update_alpha_hist")


#' @rdname update_alpha
#
#' @method update_alpha_hist multfsusie
#
#' @export update_alpha_hist.multfsusie
#
#' @export
#' @keywords internal
update_alpha_hist.multfsusie <-  function(multfsusie.obj , discard=FALSE, ... )
{
  if(!discard){
    multfsusie.obj$alpha_hist[[ (length(multfsusie.obj$alpha_hist)+1)  ]] <- multfsusie.obj$alpha
  }
  if(discard){
    if((length(multfsusie.obj$alpha_hist[[length(multfsusie.obj$alpha_hist)]]) >multfsusie.obj$L)){

      tt <- multfsusie.obj$alpha_hist[[ (length(multfsusie.obj$alpha_hist) ) ]][1:multfsusie.obj$L]
      multfsusie.obj$alpha_hist[[ (length(multfsusie.obj$alpha_hist))  ]] <- tt
    }
  }

  return( multfsusie.obj)
}



#@title Update  multfsusie object using the output of EM_pi
#
# @param multfsusie.obj a multfsusie object defined by \code{\link{init_multfsusie_obj }} function
#
# @param l integer larger or equal to 1. Corresponds to the effect to be accessed
#
# @param EM_pi an object of the class "EM_pi" generated by the function \code\link{EM_pi_multfsusie}}
# @param effect_estimate a list of marginal association generated by \code{cal_Bhat_Shat_multsusie}
# @param list_indx_lst list of list generated by \code{\link{gen_wavelet_indx}} for the given level of resolution
# @param low_trait object to check if data have critically low variance
#
# @return multfsusie object
#
# @export

update_multfsusie   <- function(multfsusie.obj, l, EM_pi, effect_estimate, list_indx_lst, low_trait,e=0.001 , ...)
{

  if( l > length(multfsusie.obj$est_pi))
  {
    stop("Error trying to access mixture proportion")
  }
  if( l < 1)
  {
    stop("Error l should be larger ")
  }
  if(  "EM_pi_multfsusie"  %!in%  class(EM_pi)  )
  {
    stop("Error EM_pi should be of the class EM_pi_multfsusie")
  }
  multfsusie.obj         <-   update_pi( multfsusie.obj =  multfsusie.obj ,
                                         l             =  l ,
                                         tpi           =  EM_pi$tpi_k)

  multfsusie.obj$G_prior <-   update_prior(get_G_prior(multfsusie.obj) , EM_pi$tpi_k  )
  post_effect <- get_post_effect_multfsusie (G_prior         = multfsusie.obj$G_prior ,
                                             effect_estimate = effect_estimate,
                                             lBF             = EM_pi$lBF,
                                             list_indx_lst   = list_indx_lst,
                                             low_trait       = low_trait,
                                             e               = e)

  if(!is.null(post_effect$res_u)){
    multfsusie.obj$ fitted_u[[l]]        <-   post_effect$res_u$Bhat
    multfsusie.obj$ fitted_u2[[l]]       <-   post_effect$res_u$Shat^2

    #fsusieR::cal_clfsr(G_prior = multfsusie.obj$G_prior$G_prior_[[k]],
    #                       Bhat = post_effect$res_u$Bhat,
    #                       Shat=post_effect$res_u$Shat,
    #                       indx_lst=list_indx_lst[[k]])

  }

  if(!is.null(post_effect$res_f)){ ### TODO: make it cleaner ------

    for (k  in 1:length(post_effect$res_f)) {

        multfsusie.obj$fitted_wc[[l]][[k]]  <- post_effect$res_f[[k]]$Bhat
        multfsusie.obj$fitted_wc2[[l]][[k]] <- post_effect$res_f[[k]]$Shat^2


    # fsusieR::cal_clfsr(G_prior = multfsusie.obj$G_prior$G_prior_f[[k]],
    #                        Bhat = post_effect$res_f[[k]]$Bhat,
    #                        Shat=post_effect$res_f[[k]]$Shat,
    #                       indx_lst=list_indx_lst[[k]])

    }
  }

  new_alpha      <- fsusieR::cal_zeta (  EM_pi$lBF)
  multfsusie.obj <- update_alpha(multfsusie.obj, l, new_alpha)
  multfsusie.obj <- update_lBF(multfsusie.obj  = multfsusie.obj,
                               l               = l,
                               lBF             = EM_pi$lBF,
                               lBF_per_trait   = EM_pi$lBF_per_trait)

 #TODO: fix that because it seems to take to much memory
 # multfsusie.obj <- update_lfsr(multfsusie.obj  = multfsusie.obj,
 #                               l               = l,
 # #                              effect_estimate = effect_estimate,
 #                               list_indx_lst   = list_indx_lst
 #                               )

  return(multfsusie.obj)
}





update_pi  <- function    (multfsusie.obj, l, tpi, ...)
  UseMethod("update_pi")

# @rdname update_pi
#
# @method update_pi multfsusie
#
# @export update_pi.multfsusie
#
#' @export
#' @keywords internal
update_pi.multfsusie <- function( multfsusie.obj, l, tpi, ...)
{

  if( l > length(multfsusie.obj$est_pi))
  {
    stop("Error trying to access mixture proportion")
  }
  if( l < 1)
  {
    stop("Error l should be larger ")
  }

  multfsusie.obj$est_pi[[l]] <- tpi
  out <- multfsusie.obj
  class(out) <- "multfsusie"
  return(out)
}


#' @title Compute KL divergence effect l
#'  @param multfsusie.obj a multfsusie object
#' @param X matrix of covariates
#
#' @param Y observed response
#
#' @param list_indx_lst  list generated by gen_wavelet_indx for the given level of resolution
#
#' @return susiF object
#' @export
#' @keywords internal


update_KL  <- function    (multfsusie.obj,Y, X , list_indx_lst , ...)
  UseMethod("update_KL")

#' @rdname update_KL
#
#' @method update_KL multfsusie
#
#' @export update_KL.multfsusie
#
#' @export
#' @keywords internal


update_KL.multfsusie <- function(multfsusie.obj, Y, X , list_indx_lst,ind_analysis, ...)
{
  multfsusie.obj$KL <-  do.call(c,lapply(1:multfsusie.obj$L,
                                         FUN=function(l)
                                           cal_KL_l(multfsusie.obj,
                                                    l=l,
                                                    Y=Y,
                                                    X=X,
                                                    list_indx_lst = list_indx_lst,
                                                    ind_analysis=ind_analysis)))
  return( multfsusie.obj)
}



#'@title Update multfsusie log Bayes factor
#
#'@param multfsusie.obj a multfsusie object defined by init_multfsusie_obj function
#'@param  ELBO new ELBO value
#'@return multfsusie object
#' @export
#' @keywords internal


update_ELBO  <- function    (multfsusie.obj,ELBO , ...)
  UseMethod("update_ELBO")


#' @rdname update_ELBO
#
#' @method update_ELBO multfsusie
#
#' @export update_ELBO.multfsusie
#
#' @export
#' @keywords internal


update_ELBO.multfsusie <- function    (multfsusie.obj,ELBO, ...)
{

  multfsusie.obj$ELBO <- c(multfsusie.obj$ELBO,ELBO)
  return(multfsusie.obj)
}




#' @title Update lfsr multfsusie l
#
#' @param multfsusie.obj a multfsusie
#'
#' @param  l the effect to be updated new ELBO value
#'
#' @return multfsusie object
#'
#' @export
#'
#' @keywords internal


update_lfsr   <- function    (multfsusie.obj, l, ...)
  UseMethod("update_lfsr")



#' @rdname update_lfsr
#'
#' @method update_lfsr multfsusie
#'
#' @export update_lfsr.multfsusie
#'
#' @export
#' @keywords internal


update_lfsr.multfsusie <- function(multfsusie.obj, l, effect_estimate, list_indx_lst,...)
{
  clfsr_mult <-  cal_clfsr(get_G_prior(multfsusie.obj),
                         effect_estimate,
                         list_indx_lst)

  if( !is.null(effect_estimate$res_f)){

      multfsusie.obj$lfsr_wc[[l]] <- lapply(1: length(effect_estimate$res_f),
                                            function(k){
                                              fsusieR:::cal_lfsr ( clfsr_mult$clfsr_wc [[k]],
                                                                       multfsusie.obj$alpha[[l]])
                                              }
                                            )



  }


  if( !is.null(effect_estimate$res_u)){

    multfsusie.obj$lfsr_u[[l]]  <- do.call( c,lapply(1: ncol(effect_estimate$res_u$Bhat),
                                                     function(k){
                                                       fsusieR:::cal_lfsr ( clfsr_mult$clfsr_u [k,],
                                                                                multfsusie.obj$alpha[[l]])
                                                     }
      )

      )



  }


  return(multfsusie.obj)
}








#' @title Update residual variance
#'
#' @description  See title
#'
#' @param multfsusie.obj a multfsusie object
#'
#' @param sigma2 the new values for residual variance
#'
#' @param ... Additional arguments passed to other functions.
#'
update_residual_variance  <- function(multfsusie.obj,sigma2, ...)
  UseMethod("update_residual_variance")



#' @rdname update_residual_variance
#
#' @method update_residual_variance multfsusie
#
#' @export update_residual_variance.multfsusie
#
#' @export
#' @keywords internal

update_residual_variance.multfsusie <- function(multfsusie.obj,sigma2, ...)
{
  multfsusie.obj$sigma2 <- sigma2
  return(multfsusie.obj)
}




#' @title Update multfsusie by computing PiP
#
#' @param multfsusie.obj a multfsusie object defined by init_multfsusie_obj function
#' @return multfsusie object
#' @export
#' @keywords internal

update_cal_pip  <- function (multfsusie.obj, ...)
  UseMethod("update_cal_pip")

#' @rdname update_cal_pip
#
#' @method update_cal_pip multfsusie
#
#' @export update_cal_pip.multfsusie
#
#' @export
#' @keywords internal

update_cal_pip.multfsusie <- function (multfsusie.obj, ...)
{
  if(sum( is.na(unlist(multfsusie.obj$alpha))))
  {
    stop("Error: some alpha value not updated, please update alpha value first")
  }
  tpip <- list()
  for ( l in 1:multfsusie.obj$L)
  {
    tpip[[l]] <- rep(1, lengths(multfsusie.obj$alpha)[[l]])-multfsusie.obj$alpha[[l]]
  }
  multfsusie.obj$pip <- 1-  apply( do.call(rbind,tpip),2, prod)
  return(multfsusie.obj)
}




#' @title Update multfsusie by computing credible sets
#
#' @param multfsusie.obj a multfsusie object defined by init_multfsusie_obj function
#
#' @param cov_lev numeric between 0 and 1, corresponding to the expected level of coverage of the cs if not specified set to 0.95
#
#' @return multfsusie object
#
#' @export
#' @keywords internal
update_cal_cs  <- function(multfsusie.obj, cov_lev=0.95, ...)
  UseMethod("update_cal_cs")

#' @rdname update_cal_cs
#
#' @method update_cal_cs multfsusie
#
#' @export update_cal_cs.multfsusie
#
#' @export
#' @keywords internal

update_cal_cs.multfsusie <- function(multfsusie.obj, cov_lev=0.95, ...)
{
  if(sum( is.na(unlist(multfsusie.obj$alpha))))
  {
    stop("Error: some alpha value not updated, please update alpha value first")
  }
  for ( l in 1:multfsusie.obj$L)
  {
    temp        <- multfsusie.obj$alpha[[l]]
    temp_cumsum <- cumsum( temp[order(temp, decreasing =TRUE)])
    max_indx_cs <- min(c(which( temp_cumsum >cov_lev ), multfsusie.obj$P))
    multfsusie.obj$cs[[l]]  <- order(temp, decreasing = TRUE)[1:max_indx_cs ]

  }

  return(multfsusie.obj)
}



#' @title Update multfsusie by computing posterior curves
#
#' @param multfsusie.obj a susiF object defined by init_multfsusie_obj function
#
#
#' @param list_indx_lst list generated by gen_wavelet_indx for the given level of resolution
#
#' @return multfsusie object
#
#' @export
#' @keywords internal
update_cal_fit_func  <- function(multfsusie.obj, list_indx_lst,...)
  UseMethod("update_cal_fit_func")

#' @rdname update_cal_fit_func
#
#' @method update_cal_fit_func multfsusie
#
#' @export update_cal_fit_func.multfsusie
#
#' @export
#' @keywords internal


update_cal_fit_func.multfsusie <- function(multfsusie.obj,list_indx_lst,... ){

  if(is.null(multfsusie.obj$fitted_wc)){
    return(multfsusie.obj)
  }

  multfsusie.obj$fitted_func <- list()
  for( l in 1: length(multfsusie.obj$cs)){

    tl <- list()

    for (k in 1: length(multfsusie.obj$n_wac))
    {
      temp <- wavethresh::wd(rep(0, multfsusie.obj$n_wac[[k]]))

      temp$D                     <- (multfsusie.obj$alpha[[l]])%*%sweep( multfsusie.obj$fitted_wc[[l]][[k]][,-list_indx_lst[[k]][[length(list_indx_lst[[k]])]]],
                                                                         1,
                                                                         1/(multfsusie.obj$csd_X ), "*")
      temp$C[length(temp$C)]     <- (multfsusie.obj$alpha[[l]])%*% (multfsusie.obj$fitted_wc[[l]][[k]][,list_indx_lst[[k]][[length(list_indx_lst[[k]])]]]*( 1/(multfsusie.obj$csd_X )))
      tl[[k]]                    <- wavethresh::wr(temp)
    }

    multfsusie.obj$fitted_func[[l]] <- tl
  }
 return( multfsusie.obj)
}



#' @title Update multfsusie by computing univariate estimates
#
#' @param multfsusie.obj a susiF object defined by init_multfsusie_obj function
#

#' @return multfsusie object
#
#' @export
#' @keywords internal
#'
update_cal_fit_u  <- function(multfsusie.obj, ...)
  UseMethod("update_cal_fit_u")

#' @rdname update_cal_fit_u
#
#' @method update_cal_fit_u multfsusie
#
#' @export update_cal_fit_u.multfsusie
#
#' @export
#' @keywords internal


update_cal_fit_u.multfsusie <- function(multfsusie.obj, ... ){

  if(is.null(multfsusie.obj$fitted_u)){
    return(multfsusie.obj)
  }


  for( l in 1: length(multfsusie.obj$cs)){

      multfsusie.obj$fitted_u[[l]] <- apply( (multfsusie.obj$alpha[[l]]) * sweep( multfsusie.obj$fitted_u[[l]] ,
                                                                         1,
                                                                         1/(multfsusie.obj$csd_X ), "*"), 2,sum)




  }
  return( multfsusie.obj)
}





#' @title Update multfsusie log Bayes factor
#
#' @param multfsusie.obj a susiF object defined by init_multfsusie_obj function
#' @param l effect to update
#' @param lBF vector of length p, containing the updated log Bayes factors
#' @return multfsusie object
#' @export
#' @keywords internal

update_lBF  <- function    (multfsusie.obj, l, lBF,...)
  UseMethod("update_lBF")


#' @rdname update_lBF
#
#' @method update_lBF multfsusie
#
#' @export update_lBF.multfsusie
#
#' @export
#' @keywords internal



update_lBF.multfsusie <- function    (multfsusie.obj,l, lBF,lBF_per_trait, ...)
{
  if(l> multfsusie.obj$L)
  {
    stop("Error: trying to update more effects than the number of specified effect")
  }

  multfsusie.obj$lBF_per_trait [[l]] <-  lBF_per_trait


  multfsusie.obj$lBF[[l]] <- lBF

  return(multfsusie.obj)
}

#' @title Check tolerance for stopping criterion
#' @description Checks whether the stopping criterion for the multfSuSiE algorithm is met.
#' @keywords internal
#' @export
test_stop_cond <- function(multfsusie.obj, check, cal_obj, Y , X, list_indx_lst, ...) {
  UseMethod("test_stop_cond")
}

#' @rdname test_stop_cond
#' @method test_stop_cond multfsusie
#' @export
test_stop_cond.multfsusie <- function(multfsusie.obj, check, cal_obj, Y, X, list_indx_lst, ind_analysis, ...) {


  if( multfsusie.obj$L==1)
  {
    multfsusie.obj$check <- 0
    return(multfsusie.obj)
  }

  if(!(multfsusie.obj$greedy_backfit_update)) #if not just updated check for stopping while loop
  {
    if( cal_obj){
      multfsusie.obj <- update_KL(multfsusie.obj,
                                  Y              = Y,
                                  X              = X,
                                  list_indx_lst  = list_indx_lst,
                                  ind_analysis   = ind_analysis)

      multfsusie.obj <- update_ELBO(multfsusie.obj,
                                    get_objective( multfsusie.obj,
                                                   Y              = Y,
                                                   X              = X,
                                                   list_indx_lst  = list_indx_lst,
                                                   ind_analysis   = ind_analysis
                                    )
      )

      if(length(multfsusie.obj$ELBO)>1    )#update parameter convergence,
      {
        check <- abs(diff(multfsusie.obj$ELBO)[(length( multfsusie.obj$ELBO )-1)])
        multfsusie.obj$check <- check
        return(multfsusie.obj)
      }else{
        multfsusie.obj$check <- check
        return(multfsusie.obj)
      }
    }
    else{
      len <- length( multfsusie.obj$alpha_hist)
      if( len>1)#update parameter convergence, no ELBO for the moment
      {
        check <-0

        T1 <- do.call( rbind, multfsusie.obj$alpha_hist[[len ]])
        T1 <- T1[1:multfsusie.obj$L,] #might be longer than L because alpha computed before discarding effect
        T2 <- do.call( rbind, multfsusie.obj$alpha_hist[[(len-1) ]])
        if(!(multfsusie.obj$L==nrow(T2))){
          return(multfsusie.obj)
        }
        T2 <- T2[1:multfsusie.obj$L,]
        if(multfsusie.obj$L==1){
          T2 <- T2[1,]
        }

        if((nrow(T1)>nrow(T2))){
          multfsusie.obj$check <- 1
          return(multfsusie.obj)
        }
        if( (nrow(T2)>nrow(T1))){
          T2 <- T2[1:multfsusie.obj$L,]
        }


        check <- sum(abs(T1-T2))/nrow(X)
        multfsusie.obj$check <- check
        return(multfsusie.obj)
        #print(check)
      }else{
        multfsusie.obj$check <- check
        return(multfsusie.obj)
      }
    }
  }else{
    multfsusie.obj$check <- check
    return(multfsusie.obj)
  }

}





#'@title Compute refined estimate using translation invariant wavelet transform
#'
#' @description   Compute refined estimate using translation invariant wavelet transform
#'
#' @param multfsusie.obj  a multfsusie.obj object
#'
#' @param Y  multfsusie.obj data object
#'
#' @param X matrix containing the covariates
#' @param verbose logical
#' @param filter.number see wd description in wavethesh package description
#' @param family  see wd description in wavethresh package description
#' @param ind_analysis output of which_notNA_pos
#'
#' @param ... Additional arguments passed to other functions.
#' @export
#' @keywords internal

TI_regression <- function (multfsusie.obj,
                           Y,
                           X,
                           ind_analysis ,
                           verbose ,
                           filter.number ,
                           family   ,...)
  UseMethod("TI_regression")


#' @rdname TI_regression
#'
#' @method TI_regression multfsusie
#'
#' @export TI_regression.multfsusie
#'
#'
#' @importFrom ashr ash
#' @importFrom fsusieR TI_regression.susiF
#' @importFrom wavethresh wd
#' @export
#' @keywords internal

TI_regression.multfsusie<- function(multfsusie.obj,
                                    Y,
                                    X,
                                    ind_analysis ,
                                    verbose=TRUE,
                                    filter.number = 10,
                                    family = "DaubLeAsymm", ...  ){

  if(is.null(multfsusie.obj$fitted_wc)){
    return(multfsusie.obj)
  }
  if(verbose){
    print( "Fine mapping done, refining effect estimates using cylce spinning wavelet transform")
  }
  tl <- list()

  for ( k in 1: length(Y$Y_f)){
    tl[[k]] <- rep (0, ncol (Y$Y_f[[k]]))


  }
  tp <- list()
  for ( l in 1:length(multfsusie.obj$cs)){
    multfsusie.obj$fitted_func[[l]] <- tl

    tp[[l]] <- list()
    for ( k in 1: length(Y$Y_f)){
      tp[[l]][[k]]  <- rbind ( rep( 0 , ncol( Y$Y_f[[k]])),
                               rep( 0 , ncol( Y$Y_f[[k]])))
    }
  }
  multfsusie.obj$cred_band <- tp
  dummy_susiF.obj <- create_dummy_susiF(multfsusie.obj)
  dummy_susiF.obj$L=multfsusie.obj$L
  for ( k in 1: length(Y$Y_f)){

    susiF.obj <- fsusieR::TI_regression.susiF(     obj           = dummy_susiF.obj,
                                                   Y             = Y$Y_f[[k]][ind_analysis$idx_f[[k]],],
                                                   X             = X[ind_analysis$idx_f[[k]],, drop=FALSE],
                                                   verbose       = FALSE,
                                                   filter.number = 1 ,
                                                   family = "DaubExPhase" )

    for ( l in 1:length(multfsusie.obj$cs)){

      multfsusie.obj$fitted_func[[l]][[k]] <- susiF.obj$fitted_func[[l]]
      multfsusie.obj$cred_band  [[l]][[k]] <- susiF.obj$cred_band  [[l]]
    }

  }

  return(multfsusie.obj)
}






#'@title Compute refined estimate using smash
#'
#' @description   Compute refined estimate using smash
#'
#' @param multfsusie.obj  a multfsusie.obj object
#'
#' @param Y  multfsusie.obj data object
#'
#' @param X matrix containing the covariates
#' @param verbose logical
#' @param filter.number see wd description in wavethesh package description
#' @param family  see wd description in wavethresh package description
#' @param ind_analysis output of which_notNA_pos
#' @param alpha required confidence level
#' @param ... Additional arguments passed to other functions.
#' @export
#' @keywords internal

smash_regression <- function (multfsusie.obj,
                           Y,
                           X,
                           ind_analysis ,
                           verbose ,
                           filter.number ,
                           family   ,...)
  UseMethod("smash_regression")


#' @rdname TI_regression
#'
#' @method TI_regression multfsusie
#'
#' @export TI_regression.multfsusie
#'
#'
#' @importFrom ashr ash
#' @importFrom fsusieR TI_regression.susiF
#' @importFrom wavethresh wd
#' @export
#' @keywords internal

smash_regression.multfsusie<- function(multfsusie.obj,
                                    Y,
                                    X,
                                    ind_analysis ,
                                    verbose=TRUE,
                                    filter.number = 10,
                                    family = "DaubLeAsymm",
                                    alpha=0.99,...   ){
  if(is.null(multfsusie.obj$fitted_wc)){
    return(multfsusie.obj)
  }
  if(verbose){
    print( "Fine mapping done, refining effect estimates using cylce spinning wavelet transform")
  }
  tl <- list()

  for ( k in 1: length(Y$Y_f)){
    tl[[k]] <- rep (0, ncol (Y$Y_f[[k]]))


  }
  tp <- list()
  for ( l in 1:length(multfsusie.obj$cs)){
    multfsusie.obj$fitted_func[[l]] <- tl

    tp[[l]] <- list()
    for ( k in 1: length(Y$Y_f)){
      tp[[l]][[k]]  <- rbind ( rep( 0 , ncol( Y$Y_f[[k]])),
                               rep( 0 , ncol( Y$Y_f[[k]])))
    }
  }
  multfsusie.obj$cred_band <- tp
  dummy_susiF.obj <- create_dummy_susiF(multfsusie.obj)

  for ( k in 1: length(Y$Y_f)){

    susiF.obj <- fsusieR::smash_regression.susiF(     obj           = dummy_susiF.obj,
                                                   Y             = Y$Y_f[[k]][ind_analysis$idx_f[[k]],],
                                                   X             = X[ind_analysis$idx_f[[k]],, drop=FALSE],
                                                   verbose       = FALSE,
                                                   filter.number = filter.number,
                                                   family        = family ,
                                                   alpha=alpha)

    for ( l in 1:length(multfsusie.obj$cs)){

      multfsusie.obj$fitted_func[[l]][[k]] <- susiF.obj$fitted_func[[l]]
      multfsusie.obj$cred_band  [[l]][[k]] <- susiF.obj$cred_band  [[l]]
    }

  }

  return(multfsusie.obj)
}





#
#' @title Return which credible sets are  dummy
#
#' @param multfsusie.obj a susif object defined by \code{\link{init_susiF_obj}} function
#' @param min_purity minimal purity within a CS
#' @param X matrix of covariate
#' @param lbf_min numeric  discard low purity cs in the IBSS fitting procedure if the largest log Bayes factors is lower than this value
#' @return a list of index corresponding the the dummy effect
#
#' @export
#' @keywords internal
#
which_dummy_cs <- function(multfsusie.obj, min_purity=0.5,X,...)
  UseMethod("which_dummy_cs")



#' @rdname which_dummy_cs
#
#' @method which_dummy_cs multfsusie
#
#' @export which_dummy_cs.multfsusie
#' @export
#' @keywords internal
#

which_dummy_cs.multfsusie  <- function(multfsusie.obj,
                                       min_purity =0.5,
                                       X,
                                       median_crit=FALSE,
                                       lbf_min,... ){
  dummy.cs<- c()
  if(  multfsusie.obj$L==1){
    return(dummy.cs)
  }
  if(missing(lbf_min)){
    lbf_min=Inf
  }
  f_crit <- function (multfsusie.obj,
                      min_purity=0.5,
                      l,
                      median_crit=FALSE,
                      lbf_min){
    if( median_crit){
      #if( length(multfsusie.obj$cs[[l]] )  < ncol(X)/10) {
      #  is.dummy.cs <- FALSE
      #   return(is.dummy.cs )
      #}
      if(length(multfsusie.obj$cs[[l]]) <5){
        is.dummy.cs <- FALSE
      }else{
        tt <-  stats::cor( X[,multfsusie.obj$cs[[l]]])

        is.dummy.cs <-   stats::median(abs( tt[lower.tri(tt, diag =FALSE)]))  <  min_purity & max(multfsusie.obj$lBF[[l]])<lbf_min
      }


    }else{
      is.dummy.cs <-   min(abs(stats::cor( X[,multfsusie.obj$cs[[l]]]))) <  min_purity & max(multfsusie.obj$lBF[[l]])<lbf_min
    }

    return( is.dummy.cs)
  }



  for (l in 1:multfsusie.obj$L )
  {

    if (length(multfsusie.obj$cs[[l]])==1)
    {

      if(   mean(unlist(get_pi0(multfsusie.obj,l=l)))>1-multfsusie.obj$tol_null_prior){# check if the estimated prior is exactly 0

        dummy.cs<-  c( dummy.cs,l)
      }

    }else{

      if(   f_crit(multfsusie.obj = multfsusie.obj,
                   min_purity=min_purity,
                   l=l,
                   median_crit=median_crit ,
                   lbf_min =lbf_min)){#check if the purity of cs l is lower that min_purity

        dummy.cs<-  c( dummy.cs,l)

      }else{
        if(   mean(unlist(get_pi0(multfsusie.obj,l=l)))>1-multfsusie.obj$tol_null_prior){

          dummy.cs<-  c( dummy.cs,l)
        }

      }
    }

  }


  if( length(dummy.cs)==0)
  {
    return(dummy.cs)
  }else{
    if(length(dummy.cs)==multfsusie.obj$L) #avoid returning empty results
    {
      dummy.cs <- dummy.cs[-length(dummy.cs)]
    }

    return(dummy.cs)
  }

}





