

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
                                         list_indx_lst  = list_indx_lst){

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


    update_D  <-  D - Reduce("+", lapply  ( id_L, function(l) (X*rep(multfsusie.obj$alpha[[l]], rep.int(dim(X)[1],dim(X)[2]))) %*% (multfsusie.obj$fitted_wc[[l]][[cord]][,-indx_lst[[length(indx_lst)]]])   ) )
    update_C  <-  C - Reduce("+", lapply  ( id_L, function(l) (X*rep(multfsusie.obj$alpha[[l]], rep.int(dim(X)[1],dim(X)[2]))) %*% multfsusie.obj$fitted_wc[[l]][[cord]][,indx_lst[[length(indx_lst)]]] ) )
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
#' @param multfsusie.obj a multfsusie object defined by \code{\link{init_multfsusie_obj}} function
#' @param min.purity minimal purity within a CS
#' @param X matrix of covariate
#' @return a multfsusie.obj without "dummy" credible s
#
#' @export
#
#
check_cs <- function(multfsusie.obj, min.purity=0.5,X,...)
  UseMethod("check_cs")



#' @rdname check_cs
#
#' @method check_cs multfsusie
#
#' @export check_cs.multfsusie
#
#' @export
#' @keywords internal


check_cs.multfsusie <- function(multfsusie.obj, min.purity=0.5,X )
{


  dummy.cs <- which_dummy_cs(multfsusie.obj, min.purity=min.purity,X)


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






#' @title Discard credible sets
#
#' @param multfsusie.obj a multfsusie object defined by \code{\link{init_multfsusie_obj}} function
#
#' @param cs vector of integer containing the credible sets to discard
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

discard_cs.multfsusie <- function(multfsusie.obj, cs, out_prep=FALSE)
{

  if( length(cs)==multfsusie.obj$L){
    cs <- cs[-1]
    if(length(cs)==0){
      return(multfsusie.obj)
    }
  }
  multfsusie.obj$alpha       <-  multfsusie.obj$alpha[ -cs]
  multfsusie.obj$lBF         <-  multfsusie.obj$lBF[ -cs]
  if(!is.null(multfsusie.obj$G_prior$G_prior_f)){
    multfsusie.obj$fitted_wc   <-  multfsusie.obj$fitted_wc[ -cs]
    multfsusie.obj$fitted_wc2  <-  multfsusie.obj$fitted_wc2[ -cs]
  }
  if(!is.null(multfsusie.obj$G_prior$G_prior_u)){
    multfsusie.obj$fitted_uni  <-  multfsusie.obj$fitted_uni[-cs]
    multfsusie.obj$fitted_uni2 <-  multfsusie.obj$fitted_uni2[-cs]
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
  L_extra <- ifelse ( max(multfsusie.obj$L_max - (multfsusie.obj$L+L_extra),0 ) >0.1,#check if we are adding more effect that maximum specified by user
                      L_extra,
                      abs(multfsusie.obj$L_max -(multfsusie.obj$L+L_extra))
                    )
  if( L_extra==0){
    return(multfsusie.obj)
  }else{
    L_old <- multfsusie.obj$L
    L_new <- multfsusie.obj$L+L_extra
    multfsusie.obj$L <- ifelse(L_new<(multfsusie.obj$P+1),L_new,P)
    l=1
    k=1

    for ( l in (L_old+1):multfsusie.obj$L )
    {

      if( !is.null(multfsusie.obj$fitted_wc)){

        multfsusie.obj$fitted_wc[[l]]  <-    lapply( 1:length(multfsusie.obj$n_wac), function(j)
                                                      matrix( 0,
                                                              ncol = multfsusie.obj$n_wac[[j]] ,
                                                              nrow =  multfsusie.obj$P))
        multfsusie.obj$fitted_wc2[[l]] <-    lapply( 1:length(multfsusie.obj$n_wac), function(j)
                                                       matrix( 1,
                                                              ncol = multfsusie.obj$n_wac[[j]],
                                                              nrow = multfsusie.obj$P))

         }
      if(!is.null(multfsusie.obj$fitted_uni)){
        multfsusie.obj$fitted_uni[[l]]        <- 0*multfsusie.obj$fitted_uni[[1]]
        multfsusie.obj$fitted_uni2[[l]]       <- 0*multfsusie.obj$fitted_uni2[[1]]

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
    return(multfsusie.obj)
  }

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
init_multfsusie_obj <- function(L_max, G_prior, Y,X,type_mark,L_start,greedy,backfit, ind_analysis )
{
  sigma2          <- list()
  if(!is.null(Y$Y_f)){
  fitted_wc       <- list()
  fitted_wc2      <- list()
  n_wac           <- lapply(lapply(Y$Y_f,dim) ,`[[`, 2)
  sigma2$sd_f     <- rep( 1, length(n_wac ))
  }else {
    fitted_wc     <- NULL
    fitted_wc2    <- NULL
    n_wac         <- NULL
    sigma2$sd_f   <- NULL
  }
  if(!is.null(Y$Y_u)){
    fitted_uni    <- list()
    fitted_uni2   <- list()

  }else{
    fitted_uni    <- NULL
    fitted_uni2   <- NULL
    sigma2$sd_u   <- NULL
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
      fitted_wc[[l]]        <-    lapply( 1:length(n_wac), function(j) matrix( 0,ncol= n_wac[[j]],  nrow = ncol(X)))
      fitted_wc2[[l]]       <-    lapply( 1:length(n_wac), function(j) matrix(  1,ncol= n_wac[[j]],  nrow = ncol(X)))
      sigma2
    }
    if(!is.null(Y$Y_u)){
      fitted_uni [[l]]       <-     matrix(0, ncol= ncol(Y$Y_u), nrow = ncol(X))
      fitted_uni2[[l]]       <-     matrix(0, ncol= ncol(Y$Y_u), nrow = ncol(X))
      sigma2$sd_u            <-   rep( 1, ncol(Y$Y_u))
    }


    alpha [[l]]           <-  rep(0, dim(X)[2])
    cs[[l]]               <-  list()
    est_pi [[l]]          <-  get_pi_G_prior(G_prior)
    lBF[[l]]              <-  rep(NA, ncol(X))

  }

 sigma2 <-  init_var_multf(Y)

  obj <- list( fitted_wc       = fitted_wc,
               fitted_wc2      = fitted_wc2,
               fitted_uni      = fitted_uni,
               fitted_uni2     = fitted_uni2,
               lBF             = lBF,
               KL              = KL,
               ELBO            = ELBO,
               ind_fitted_val  = ind_fitted_val,
               G_prior         = G_prior,
               alpha_hist      = alpha_hist,
               N               = N,
               n_wac           = n_wac,
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
               ind_analysis    =  ind_analysis)

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

        pi0_f <-lapply( 1:length(multfsusie.obj$n_wac),
                   function(k)
                    sapply(multfsusie.obj$est_pi[[l]]$est_pi_f[[k]],"[[",1)
          )

        out[[l]]  <- list( pi0_u=pi0_u,
                           pi0_f = pi0_f)
    }



    }
  }else{

    pi0_u <- lapply( 1:length(multfsusie.obj$n_wac),
                     function(k)
                       sapply(multfsusie.obj$est_pi[[l]]$est_pi_u[[k]],"[[",1)
              )
    pi0_f <- lapply( 1:length(multfsusie.obj$n_wac),
                     function(k)
                       sapply(multfsusie.obj$est_pi[[l]]$est_pi_f[[k]],"[[",1)
             )


      out   <- list( pi0_u=pi0_u,
                         pi0_f = pi0_f)
    }

  return(out)
}




#' @title Acces prior object for multfsusie.obj
#' @description see title
#' @param multfsusie.obj a multfsusie object
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
get_G_prior.multfsusie <- function(multfsusie.obj){
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
get_lBF.multfsusie <- function(multfsusie.obj,l){
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
get_post_effect_multfsusie <- function(G_prior, effect_estimate, list_indx_lst=NULL, low_trait = NULL){

  out <- list( res_uni = NULL,
               res_f  = NULL)
 if( !is.null(effect_estimate$res_uni)){


     out$res_uni$Bhat  <- do.call(
                                  cbind,
                                  lapply( 1:length(G_prior$G_prior_u),
                                          function( k) get_post_mean_u(G_prior$G_prior_u[[k]],
                                                                       effect_estimate$res_uni$Bhat[,k],
                                                                       effect_estimate$res_uni$Shat[,k],
                                                                       low_u = ifelse(k%in%low_trait$low_u, TRUE,FALSE)
                                                                        )
                                          )
                                  )
     out$res_uni$Shat  <- do.call(
                                  cbind,
                                  lapply( 1:length(G_prior$G_prior_u),
                                          function( k) get_post_sd_u(G_prior$G_prior_u[[k]],
                                                                     effect_estimate$res_uni$Bhat[,k],
                                                                     effect_estimate$res_uni$Shat[,k],
                                                                     low_u = ifelse(k%in%low_trait$low_u, TRUE,FALSE)
                                                                        )
                                         )
                                  )

   }

  if( !is.null(effect_estimate$res_f)){
   out$res_f <- lapply( 1:length(G_prior$G_prior_f), function(k)list_post_mean_sd(G_prior$G_prior_f[[k]],
                                                               effect_estimate$res_f[[k]]$Bhat,
                                                               effect_estimate$res_f[[k]]$Shat,
                                                               list_indx_lst[[k]],
                                                               lowc_wc= low_trait$low_wc[[k]] )
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
# @export
#
get_post_F.multfsusie <- function (multfsusie.obj,l){

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
    if(!is.null(multfsusie.obj$fitted_uni)){
      out_u <-  Reduce("+",
                       lapply(1:multfsusie.obj$L,
                              function(l)
                                multfsusie.obj$alpha[[l]] * multfsusie.obj$fitted_uni[[l]]
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
      if(!is.null(multfsusie.obj$fitted_uni)){
             out_u <-  multfsusie.obj$alpha[[l]] * multfsusie.obj$fitted_uni[[l]]
      }else{
             out_u <- NULL
      }
  }

  out_post <- list(post_uni = out_u,
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
#
get_post_F2.multfsusie <- function (multfsusie.obj,l){

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
    if(!is.null(multfsusie.obj$fitted_uni)){
      out_u <-  Reduce("+",
                       lapply(1:multfsusie.obj$L,
                              function(l)
                                multfsusie.obj$alpha[[l]] * (multfsusie.obj$fitted_uni2[[l]]+multfsusie.obj$fitted_uni[[l]]^2    )
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
    if(!is.null(multfsusie.obj$fitted_uni)){
      out_u <-   multfsusie.obj$alpha[[l]] * (multfsusie.obj$fitted_uni2[[l]]+multfsusie.obj$fitted_uni[[l]]^2    )

    }else{
      out_u <- NULL
    }
  }

  out_post <- list(post_uni_sd2 = out_u,
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

get_ER2.multfsusie = function (  multfsusie.obj,Y, X,ind_analysis ) {
  postF <- get_post_F(multfsusie.obj )# J by N matrix
  #Xr_L = t(X%*% postF)
  postF2 <- get_post_F2(multfsusie.obj ) # Posterior second moment.

  ER2 <-  list()
  if(! is.null(Y$Y_u)){

    if( missing(ind_analysis)){
      ER2$uni <-  do.call( c,
                           lapply(1:ncol( Y$Y_u),
                                  function(k) sum((Y$Y_u[,k] - X%*%postF$post_uni[,k] )^2) -sum( postF$post_uni_sd2[,k]^2) +sum( postF2$post_uni_sd2[,k])
                           )
      )
    }else{
      ER2$uni <-  do.call( c,
                           lapply(1:ncol( Y$Y_u),
                                  function(k) sum((Y$Y_u[ind_analysis$idx_u[[k]],k] - X[ind_analysis$idx_u[[k]],]%*%postF$post_uni[,k] )^2) -sum( postF$post_uni_sd2[,k]^2) +sum( postF2$post_uni_sd2[,k])
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
                                function(k) sum((Y$Y_f[[k]][ind_analysis$idx_f[[k]],] - X[ind_analysis$idx_f[[k]],]%*%postF$post_f[[k]])^2)  -sum(postF$post_f[[k]]^2) + sum( postF2$post_f_sd2 [[k]])
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
#' @param multfsusie.obj a susiF object defined by \code{\link{init_multfsusie.obj_obj}} function
#
#' @return multfsusie.obj object
#
#' @importFrom susiF.alpha cal_cor_cs
#' @export get_pi0.multfsusie
#' @keywords internal
#
greedy_backfit  <-  function(multfsusie.obj, verbose,cov_lev,X,min.purity, ...  )
  UseMethod("greedy_backfit")

#' @rdname greedy_backfit
#'
#' @method greedy_backfit multfsusie
#'
#' @export greedy_backfit.multfsusie
#
#' @export get_pi0.multfsusie
#' @keywords internal

greedy_backfit.multfsusie <-  function(multfsusie.obj,verbose,cov_lev,X,min.purity, ...  )
{


  multfsusie.obj <- update_alpha_hist(multfsusie.obj)
  if(!(multfsusie.obj$greedy)&!(multfsusie.obj$backfit))
  {
    return(multfsusie.obj)
  }
  multfsusie.obj <- update_cal_cs(multfsusie.obj,
                                  cov_lev=cov_lev)

  dummy.cs <-  which_dummy_cs(multfsusie.obj,
                              min.purity = min.purity,
                              X=X)
  if(multfsusie.obj$backfit & (length(dummy.cs)>0)){

    multfsusie.obj$greedy <- FALSE
    if(length(dummy.cs)== multfsusie.obj$L){
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

      if( length(multfsusie.obj$cs)>1){
        A <- susiF.alpha::cal_cor_cs(multfsusie.obj, X)$cs_cor
        tl <- which(A>0.99, arr.ind = TRUE)
        tl <-  tl[- which( tl[,1]==tl[,2]),]

        if ( length(tl )==0){

        }else{

          tl <-  tl[which(tl[,1] < tl[,2]),]
          multfsusie.obj <- merge_effect(multfsusie.obj, tl)

        }
      }
      if(verbose){
        print( paste( "Discarding ",(temp_L- multfsusie.obj$L), " effects"))
      }
    }
    return(multfsusie.obj)

  }##Conditions for stopping greedy search
  if(  (multfsusie.obj$L>multfsusie.obj$L_max))
  {

    multfsusie.obj$greedy <- FALSE



    multfsusie.obj <- discard_cs(multfsusie.obj,
                                 cs= (multfsusie.obj$L_max+1):multfsusie.obj$L,
                                 out_prep= FALSE
    )




    if( length(multfsusie.obj$cs)>1){
      A <- susiF.alpha::cal_cor_cs(multfsusie.obj, X)$cs_cor
      tl <- which(A>0.99, arr.ind = TRUE)
      tl <-  tl[- which( tl[,1]==tl[,2]),]

      if ( dim(tl)[1]==0){

      }else{

        tl <-  tl[which(tl[,1] < tl[,2]),]
        multfsusie.obj <- merge_effect(multfsusie.obj, tl)

      }
    }
    if(verbose){
      print( paste( "Discarding ",(multfsusie.obj$L_max- multfsusie.obj$L), " effects"))
      print( "Greedy search and backfitting done")
    }

  }

  if( length(dummy.cs)==0& !( multfsusie.obj$greedy))
  {

    multfsusie.obj$backfit <- FALSE
  }

  if(!(multfsusie.obj$greedy )&!(multfsusie.obj$backfit ) ){

    if( length(multfsusie.obj$cs)>1){
      A <- susiF.alpha::cal_cor_cs(multfsusie.obj, X)$cs_cor
      tl <- which(A>0.99, arr.ind = TRUE)
      tl <-  tl[- which( tl[,1]==tl[,2]),]

      if ( dim(tl)[1]==0){

      }else{

        tl <-  tl[which(tl[,1] < tl[,2]),]
        multfsusie.obj <- merge_effect(multfsusie.obj, tl)

      }
    }

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
      if( length(multfsusie.obj$cs)>1){
        A <- susiF.alpha::cal_cor_cs(multfsusie.obj, X)$cs_cor
        tl <- which(A>0.99, arr.ind = TRUE)
        tl <-  tl[- which( tl[,1]==tl[,2]),]

        if ( dim(tl)[1]==0){

        }else{

          tl <-  tl[which(tl[,1] < tl[,2]),]
          multfsusie.obj <- merge_effect(multfsusie.obj, tl)

        }
      }

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



    if( length(multfsusie.obj$cs)>1){
      A <- susiF.alpha::cal_cor_cs(multfsusie.obj, X)$cs_cor
      tl <- which(A>0.99, arr.ind = TRUE)
      tl <-  tl[- which( tl[,1]==tl[,2]),]

      if ( dim(tl)[1]==0){

      }else{

        tl <-  tl[which(tl[,1] < tl[,2]),]
        multfsusie.obj <- merge_effect(multfsusie.obj, tl, discard=FALSE)

      }
    }



    multfsusie.obj <- expand_multfsusie_obj(multfsusie.obj,L_extra = temp)
    return(multfsusie.obj)
  }

}



# @title formatting function for output of posterior quantities
# @description formatting function for output of posterior quantities
# @param G_prior a mixutre_per_scale prior
# @param  Bhat matrix of estimated mean
# @param Shat matrix of estimated sd
# @importFrom susiF.alpha post_mat_mean
# @importFrom susiF.alpha post_mat_sd
#

list_post_mean_sd <- function(G_prior, Bhat,Shat,  indx_lst, lowc_wc=NULL)
{
  out <- list (Bhat= susiF.alpha::post_mat_mean( G_prior ,
                                                 Bhat,
                                                 Shat,
                                                 indx_lst,
                                                 lowc_wc=lowc_wc),
               Shat=susiF.alpha:: post_mat_sd(   G_prior ,
                                                 Bhat,
                                                 Shat,
                                                 indx_lst,
                                                 lowc_wc=lowc_wc)
  )
  return(out)

}



#' @title Merging effect function
#
#' @param multfsusie.obj a susiF object defined by \code{\link{init_multfsusie_obj}} function
#
#' @param tl see  \code{\link{greedy_backfit}}
#
#
#
#' @return  a multfusie  object
#' @export
#' @keywords internal
merge_effect <- function( multfsusie.obj, tl, ...)
  UseMethod("merge_effect")

#' @rdname merge_effect
#
#' @method merge_effect multfsusie
#
#' @export merge_effect.multfsusie
#
#' @export
#' @keywords internal

merge_effect.multfsusie <- function( multfsusie.obj, tl, discard=TRUE){




  if(is.vector( tl)){
    #print( tl)
    if( !is.null(multfsusie.obj$fitted_wc[[1]])){
      for ( k in 1: length(multfsusie.obj$fitted_wc[[1]])){
        multfsusie.obj$fitted_wc[[tl[  2]]][[k]] <- 0* multfsusie.obj$fitted_wc[[tl[ 2]]][[k]]
        multfsusie.obj$fitted_wc[[tl[  1]]][[k]] <- multfsusie.obj$fitted_wc[[tl[  1]]][[k]]+   multfsusie.obj$fitted_wc[[tl[ 2]]][[k]]
        multfsusie.obj$fitted_wc2[[tl[ 1]]][[k]] <- multfsusie.obj$fitted_wc2[[tl[  1]]][[k]] +   multfsusie.obj$fitted_wc2[[tl[  2]]][[k]]

      }


    }
    if(!is.null(multfsusie.obj$fitted_uni[[1]])){

        multfsusie.obj$fitted_uni [[tl[  2]]]  <- 0* multfsusie.obj$fitted_uni[[tl[ 2]]]
        multfsusie.obj$fitted_uni [[tl[  1]]]  <- multfsusie.obj$fitted_uni[[tl[  1]]]  +   multfsusie.obj$fitted_uni[[tl[ 2]]]
        multfsusie.obj$fitted_uni2[[tl[  1]]] <- multfsusie.obj$fitted_uni2[[tl[  1]]]  +   multfsusie.obj$fitted_uni2[[tl[  2]]]

    }
    tindx <-  tl[  2]
  }else{
    tl <- tl[order(tl[,1], tl[,2], decreasing = TRUE),]
    #print( tl)
    tindx <- c(0)
    for ( o in 1:dim(tl)[1]){

      if ( tl[o, 2]%!in%tindx){
        if( !is.null(multfsusie.obj$fitted_wc[[1]])){
          for ( k in 1: length(multfsusie.obj$fitted_wc[[1]])){
            multfsusie.obj$fitted_wc[[tl[  2]]][[k]] <- 0* multfsusie.obj$fitted_wc[[tl[ 2]]][[k]]
            multfsusie.obj$fitted_wc[[tl[  1]]][[k]] <- multfsusie.obj$fitted_wc[[tl[  1]]][[k]] +   multfsusie.obj$fitted_wc[[tl[ 2]]][[k]]
            multfsusie.obj$fitted_wc2[[tl[  1]]][[k]] <- multfsusie.obj$fitted_wc2[[tl[  1]]][[k]] +   multfsusie.obj$fitted_wc2[[tl[  2]]][[k]]

          }


        }
        if(!is.null(multfsusie.obj$fitted_uni[[1]])){

          multfsusie.obj$fitted_uni[[tl[   2]]]  <- 0* multfsusie.obj$fitted_uni[[tl[ 2]]]
          multfsusie.obj$fitted_uni[[tl[   1]]]  <- multfsusie.obj$fitted_uni[[tl[    1]]]  +   multfsusie.obj$fitted_uni[[tl[ 2]]]
          multfsusie.obj$fitted_uni2[[tl[  1]]] <- multfsusie.obj$fitted_uni2[[tl[    1]]]  +   multfsusie.obj$fitted_uni2[[tl[  2]]]

        }
        tindx <- c(tindx, tl[o, 2])
      }

    }

    tindx <- tindx[-1]

  }
  if(discard){
    multfsusie.obj<-  discard_cs(multfsusie.obj,cs=tindx, out_prep=FALSE)
  }

  return( multfsusie.obj)
}



# @title Updates CS names for output
#
# @param multfsusie.obj a multfsusie object defined by \code{\link{init_multfsusie_obj}} function
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
# @export
#

name_cs.multfsusie <- function(multfsusie.obj,X,...){

  if( length(colnames(X))==ncol(X)){

    for (l in 1: length(multfsusie.obj$cs)){
      names(multfsusie.obj$cs[[l]]) <- colnames(X)[multfsusie.obj$cs[[l]]]
    }

  }
  return(multfsusie.obj)
}






pred_partial_u <- function( multfsusie.obj, l, X )
{
  tbeta <-multfsusie.obj$alpha[[l]] *multfsusie.obj$fitted_uni[[l]]
  pred_l   <- X%*% ( tbeta)
  return(pred_l)
}

#' @title Update alpha   susiF mixture proportion of effect l
#
#' @param multfsusie.obj a multfsusie  object defined by \code{\link{init_multfsusie_obj}} function
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
#' @param multfsusie.obj a multfsusie object defined by \code{\link{init_multfsusie_obj}} function
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

update_multfsusie   <- function(multfsusie.obj, l, EM_pi, effect_estimate, list_indx_lst, low_trait,  ...)
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
                                             list_indx_lst   = list_indx_lst,
                                             low_trait       = low_trait)

  if(!is.null(post_effect$res_uni)){
    multfsusie.obj$ fitted_uni[[l]]        <-   post_effect$res_u$Bhat
    multfsusie.obj$ fitted_uni2[[l]]       <-   post_effect$res_u$Shat^2
  }

  if(!is.null(post_effect$res_f)){ ### TODO: make it cleaner ------

    for (k  in 1:length(post_effect$res_f)) {

      if( is.null(low_trait$low_wc[[k]])){
        multfsusie.obj$fitted_wc[[l]][[k]]  <- post_effect$res_f[[k]]$Bhat
        multfsusie.obj$fitted_wc2[[l]][[k]] <- post_effect$res_f[[k]]$Shat^2
      }else{
        multfsusie.obj$fitted_wc[[l]][[k]][,-low_trait$low_wc[[k]]]    <-  post_effect$res_f[[k]]$Bhat
        multfsusie.obj$fitted_wc2[[l]][[k]][,-low_trait$low_wc[[k]]]   <- post_effect$res_f[[k]]$Shat^2
      }
    }
  }

  new_alpha      <- susiF.alpha::cal_zeta (  EM_pi$lBF)
  multfsusie.obj <- update_alpha(multfsusie.obj, l, new_alpha)
  multfsusie.obj <- update_lBF(multfsusie.obj, l, EM_pi$lBF)

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
# @export
#
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
                                           cal_KL_l(multfsusie.obj, l, Y, X,
                                                    list_indx_lst = list_indx_lst,
                                                    ind_analysis=ind_analysis)))
  return( multfsusie.obj)
}



#'@title Update multfsusie log Bayes factor
#
#'@param multfsusie.obj a multfsusie object defined by \code{\link{init_multfsusie_obj}} function
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






#' @title Update residual variance
#' @description  See title
#' @param multfsusie.obj a multfsusie object
#' @param sigma2 the new values for residual variance
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

update_residual_variance.multfsusie <- function(multfsusie.obj,sigma2)
{
  multfsusie.obj$sigma2 <- sigma2
  return(multfsusie.obj)
}




#' @title Update multfsusie by computing PiP
#
#' @param multfsusie.obj a multfsusie object defined by \code{\link{init_multfsusie_obj}} function
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
#' @param multfsusie.obj a multfsusie object defined by \code{\link{init_multfsusie_obj}} function
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

update_cal_cs.multfsusie <- function(multfsusie.obj, cov_lev=0.95)
{
  if(sum( is.na(unlist(multfsusie.obj$alpha))))
  {
    stop("Error: some alpha value not updated, please update alpha value first")
  }
  for ( l in 1:multfsusie.obj$L)
  {
    temp        <- multfsusie.obj$alpha[[l]]
    temp_cumsum <- cumsum( temp[order(temp, decreasing =TRUE)])
    max_indx_cs <- min(which( temp_cumsum >cov_lev ))
    multfsusie.obj$cs[[l]]  <- order(temp, decreasing = TRUE)[1:max_indx_cs ]

  }

  return(multfsusie.obj)
}


#' @title Preparing output of main multfsusie function
#
#' @param multfsusie.obj a multfsusie object defined by \code{\link{init_multfsusie_obj}} function
#
#' @param Y  data
#' @param X matrix of size N by p
#
#' @param list_indx_lst list generated by gen_wavelet_indx for the given level of resolution
#
#' @param filter.cs logical, if TRUE filter the credible set (removing low purity cs and cs with estimated prior equal to 0)
#' @export
#' @keywords internal
out_prep <- function(multfsusie.obj,Y,X,list_indx_lst,filter.cs, ...)
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
                                X,
                                list_indx_lst,
                                filter.cs,
                                outing_grid)
{
  multfsusie.obj <-  update_cal_pip(multfsusie.obj)
  multfsusie.obj <-  update_cal_fit_func(multfsusie.obj,list_indx_lst)

  multfsusie.obj <-  update_cal_fit_uni(multfsusie.obj )

  multfsusie.obj <-  name_cs(multfsusie.obj,X)
  if(filter.cs)
  {
    multfsusie.obj <- check_cs(multfsusie.obj,min.purity=0.5,X=X)
  }
  multfsusie.obj$outing_grid <-outing_grid
  multfsusie.obj$purity      <-susiF.alpha::cal_purity(l_cs= multfsusie.obj$cs, X=X)
  return( multfsusie.obj)
}


#' @title Update multfsusie by computing posterior curves
#
#' @param multfsusie.obj a susiF object defined by \code{\link{init_multfsusie_obj}} function
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
#' @param multfsusie.obj a susiF object defined by \code{\link{init_multfsusie_obj}} function
#

#' @return multfsusie object
#
#' @export
#' @keywords internal
#'
update_cal_fit_uni  <- function(multfsusie.obj, ...)
  UseMethod("update_cal_fit_uni")

#' @rdname update_cal_fit_uni
#
#' @method update_cal_fit_uni multfsusie
#
#' @export update_cal_fit_uni.multfsusie
#
#' @export
#' @keywords internal


update_cal_fit_uni.multfsusie <- function(multfsusie.obj, ... ){

  if(is.null(multfsusie.obj$fitted_uni)){
    return(multfsusie.obj)
  }


  for( l in 1: length(multfsusie.obj$cs)){

      multfsusie.obj$fitted_uni[[l]] <- (multfsusie.obj$alpha[[l]]) * sweep( multfsusie.obj$fitted_uni[[l]] ,
                                                                         1,
                                                                         1/(multfsusie.obj$csd_X ), "*")




  }
  return( multfsusie.obj)
}





#' @title Update multfsusie log Bayes factor
#
#' @param multfsusie.obj a susiF object defined by \code{\link{init_multfsusie_obj}} function
#' @param l effect to update
#' @param lBF vector of length p, containing the updated log Bayes factors
#' @return multfsusie object
#' @export
#' @keywords internal

update_lBF  <- function    (susiF.obj, l, lBF,...)
  UseMethod("update_lBF")


#' @rdname update_lBF
#
#' @method update_lBF multfsusie
#
#' @export update_lBF.multfsusie
#
#' @export
#' @keywords internal



update_lBF.multfsusie <- function    (multfsusie.obj,l, lBF, ...)
{
  if(l> multfsusie.obj$L)
  {
    stop("Error: trying to update more effects than the number of specified effect")
  }

  multfsusie.obj$lBF[[l]] <- lBF
  return(multfsusie.obj)
}

#
#' @title Check tolerance for stopping criterion
#
#' @export
#' @keywords internal
#
test_stop_cond <- function(multfsusie.obj, check, cal_obj, Y_data, X, list_indx_lst  ,...)
  UseMethod("test_stop_cond")



#' @rdname test_stop_cond
#
#' @method test_stop_cond multfsusie
#
#' @export test_stop_cond.multfsusie
#
#' @export
#


test_stop_cond.multfsusie<- function(multfsusie.obj, check, cal_obj, Y, X, list_indx_lst, ind_analysis  ,...)
{

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




#
#' @title Return which credible sets are  dummy
#
#' @param multfsusie.obj a susif object defined by \code{\link{init_susiF_obj}} function
#' @param min.purity minimal purity within a CS
#' @param X matrix of covariate
#
#' @return a list of index corresponding the the dummy effect
#
#' @export
#' @keywords internal
#
which_dummy_cs <- function(multfsusie.obj, min.purity=0.5,X,...)
  UseMethod("which_dummy_cs")



#' @rdname which_dummy_cs
#
#' @method which_dummy_cs multfsusie
#
#' @export which_dummy_cs.multfsusie
#' @export
#' @keywords internal
#

which_dummy_cs.multfsusie  <- function(multfsusie.obj, min.purity =0.5, X){
  dummy.cs<- c()
  if(  multfsusie.obj$L==1){
    return(dummy.cs)
  }


  for (l in 1:multfsusie.obj$L )
  {

    if (length(multfsusie.obj$cs[[l]])==1)
    {

      if(   mean(unlist(get_pi0(multfsusie.obj,l=l)))==1){# check if the estimated prior is exactly 0

        dummy.cs<-  c( dummy.cs,l)
      }

    }else{

      if( min(cor( X[,multfsusie.obj$cs[[l]]])) <  min.purity){#check if the purity of cs l is lower that min.purity

        dummy.cs<-  c( dummy.cs,l)

      }else{
        if(  mean(unlist(get_pi0(multfsusie.obj,l=l)))==1){

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
