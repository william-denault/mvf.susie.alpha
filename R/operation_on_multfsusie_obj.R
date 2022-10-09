
#' @rdname cal_partial_resid
#'
#' @method cal_partial_resid multfsusie
#'
#' @export cal_partial_resid.multfsusie
#'
#' @export
#'


cal_partial_resid.multfsusie <- function(multfsusie.obj = multfsusie.obj,l = l, X = X,Y = Y,list_indx_lst  = list_indx_lst){

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
    update_Y$Y_f   <-   lapply(1:length(Y$Y_f),function(k)   cal_partial_resid_sub(multfsusie.obj,
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



#'
#' @title Check purity credible sets
#'
#' @param multfsusie.obj a multfsusie object defined by \code{\link{init_multfsusie_obj}} function
#' @param min.purity minimal purity within a CS
#' @param X matrix of covariate
#' @return a multfsusie.obj without "dummy" credible s
#'
#' @export
#'
#'
check_cs <- function(multfsusie.obj, min.purity=0.5,X,...)
  UseMethod("check_cs")



#' @rdname check_cs
#'
#' @method check_cs multfsusie
#'
#' @export check_cs.multfsusie
#'
#' @export
#'

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
#'
#' @param multfsusie.obj a multfsusie object defined by \code{\link{init_multfsusie_obj}} function
#'
#' @param cs vector of integer containing the credible sets to discard
#'
#' @return a multfsusie.obj without "dummy" credible sets
#'
#' @export



discard_cs <- function(multfsusie.obj, cs,...)
  UseMethod("discard_cs")



#' @rdname discard_cs
#'
#' @method discard_cs multfsusie
#'
#' @export discard_cs.multfsusie
#'
#' @export
#'

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
    multfsusie.obj$ELBO                  <- c()
  }

  multfsusie.obj$cs          <-  multfsusie.obj$cs[ -cs]
  multfsusie.obj$est_pi           <-  multfsusie.obj$est_pi[ -cs]
  multfsusie.obj$L                <-  multfsusie.obj$L -length(cs)
  multfsusie.obj$est_pi      <- multfsusie.obj$est_pi[ -cs]
  return(multfsusie.obj)
}




#' @title Expand multfsusie.obj by adding L_extra effect
#'
#' @param multfsusie.obj a multfsusie.obj
#'
#' @param L_extra numeric a number of effect to add
#'
#' @return a multfsusie.obj a L_extra effect. Note the the number of effect of the multfsusie.obj cannot exceed the number the user upper bound
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

    for ( l in (L_old+1):multfsusie.obj$L )
    {
      multfsusie.obj$fitted_wc[[l]]        <-  lapply(1:length(multfsusie.obj$fitted_wc[[1]]),
                                                      function (k) 0*multfsusie.obj$fitted_wc[[1]][k])
      multfsusie.obj$fitted_wc2[[l]]       <-   lapply(1:length(multfsusie.obj$fitted_wc[[1]]),
                                                       function (k) (0*multfsusie.obj$fitted_wc2[[1]][k] +1))
      multfsusie.obj$alpha [[l]]           <-  rep(0, length(multfsusie.obj$alpha [[1]]))
      multfsusie.obj$cs[[l]]               <-  list()
      multfsusie.obj$est_pi [[l]]          <-  multfsusie.obj$est_pi[[1]]
      #multfsusie.obj$est_sd [[l]]          <-  multfsusie.obj$est_sd[[1]]
      multfsusie.obj$lBF[[l]]              <-  rep(NA, length( multfsusie.obj$lBF[[1]]))
      multfsusie.obj$KL                    <- rep(NA,multfsusie.obj$L)
      multfsusie.obj$ELBO                  <- c()
    }
    multfsusie.obj$n_expand <- multfsusie.obj$n_expand+1
    multfsusie.obj$greedy_backfit_update <- TRUE
    return(multfsusie.obj)
  }

}



#' @title Initialize a multsusie object
#'
#' @param L_max upper bound on the number of non zero coefficients An L-vector containing the indices of the
#'   nonzero coefficients.
#'
#'
#' @param G_prior prior object defined by init_prior_multsusie function
#'
#' @param Y list of matrices of outcomes
#'
#' @param X Matrix of covariates
#'
#' @param L_start number of effect to start with
#' @pama type_mark an object generated by \code{\link{is.functional}} function
#'
#' @return A list with the following elements
#' \item{fitted_wc}{ list of length L, each element contains the fitted wavelet coefficients of effect l}
#' \item{fitted_wc2}{list of length L, each element contains the variance of the fitted wavelet coefficients of effect l}
#' \item{alpha_hist}{ history of the fitted alpha value}
#' \item{N}{ number of indidivual in the study}
#' \item{sigma2}{residual variance}
#' \item{n_wac}{number of wavelet coefficients}
#' \item{ind_fitted_func}{fitted curves of each individual }
#' \item{cs}{credible set}
#' \item{pip}{Posterior inclusion probabilites}
#' \item{G_prior}{a G_prior of the same class as the input G_prior, used for internal calculation}
#' \item{lBF}{ log Bayes factor for the different effect}
#' \item{KL}{ the KL divergence for the different effect}
#' \item{ELBO}{ The evidence lower bound}
#' \item{lfsr_wc}{Local fasle sign rate of the fitted wavelet coefficients}
#' @export
init_multfsusie_obj <- function(L_max, G_prior, Y,X, type_mark,L_start,greedy,backfit )
{
  sigma2          <-  list()
  if(!is.null(Y$Y_f)){
  fitted_wc       <-  list()
  fitted_wc2      <-  list()
  n_wac           <-  lapply(lapply(Y$Y_f,dim) ,`[[`, 2)
  sigma2$sd_f     <-  rep( 1, length(n_wac ))
  }else {
    fitted_wc       <-  NULL
    fitted_wc2      <-  NULL
    n_wac           <-  NULL
    sigma2$sd_f     <-  NULL
  }
  if(!is.null(Y$Y_u)){
    fitted_uni        <-   list()
    fitted_uni2       <-   list()

  }else{
    fitted_uni        <-   NULL
    fitted_uni2       <-   NULL
    sigma2$sd_u       <-   NULL
  }
  alpha           <-  list()
  alpha_hist      <-  list()
  ind_fitted_val  <-  list()
  cs              <-  list()
  pip             <-  rep(0, dim(X)[2])
  est_pi          <-  list()
  est_sd          <-  list()
  L_max           <-  L_max
  L               <-  min(3,L_max)
  G_prior         <-  G_prior
  N               <-  nrow(X)[1]
  n_cond          <-  type_mark$ncond
  P               <-  ncol(X)
  lBF             <-  list()
  KL              <-  rep(NA,L)
  ELBO            <-  c()
  mean_X          <- attr(X, "scaled:center")
  csd_X           <- attr(X, "scaled:scale")
  n_expand        <- 0 #number of greedy expansion
  greedy          <- greedy
  backfit         <- backfit
  greedy_backfit_update <- FALSE
  for ( l in 1:L )
  {

    if(!is.null(Y$Y_f)){
      fitted_wc[[l]]        <-    lapply( 1:length(n_wac), function(j) matrix( 0,ncol= n_wac[[j]],  nrow = ncol(X)))
      fitted_wc2[[l]]       <-    lapply( 1:length(n_wac), function(j) rep(  0,ncol= n_wac[[j]],  nrow = ncol(X)))
      sigma2
    }
    if(!is.null(Y$Y_u)){
      fitted_uni [[l]]       <-    matrix(0, ncol= ncol(Y$Y_u), nrow = ncol(X))
      fitted_uni2[[l]]       <-     matrix(0, ncol= ncol(Y$Y_u), nrow = ncol(X))
      sigma2$sd_u            <-   rep( 1, ncol(Y$Y_u))
    }


    alpha [[l]]           <-  rep(0, dim(X)[2])
    cs[[l]]               <-  list()
    est_pi [[l]]          <-  get_pi_G_prior(G_prior)
    lBF[[l]]              <-  rep(NA, ncol(X))

  }



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
               L               = L,
               L_max           = L_max,
               greedy          = greedy,
               backfit         = backfit,
               greedy_backfit_update=greedy_backfit_update)

  class(obj) <- "multfsusie"
  return(obj)
}







get_pi <- function(multfsusie.obj,l, ...)
  UseMethod("get_pi")



#' @rdname get_pi
#'
#' @method get_pi multfsusie
#'
#' @export get_pi.multfsusie
#'
#' @export
#'
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




#' @rdname get_G_prior
#'
#' @method get_G_prior multfsusie
#'
#' @export get_G_prior.multfsusie
#'
#' @export
#'
get_G_prior.multfsusie <- function(multfsusie.obj){
  out <- multfsusie.obj$G_prior
  return(out)
}



#' @rdname get_lBF
#'
#' @method get_lBF multfsusie
#'
#' @export get_lBF.multfsusie
#'
#' @export
#'
get_lBF.multfsusie <- function(multfsusie.obj,l){
  out <- multfsusie.obj$lBF[[l]]
  return(out)
}



#' @title Compute posterior mean of the fitted effect
#'
#' @param G_prior a multfsusie_prior object
#'
#' @param effect_estimate an object generated by \link{\code{cal_Bhat_Shat_multfsusie }}
#' @param list_indx_lst list of list generated by \code{\link{gen_wavelet_indx}} for the given level of resolution
#' @param lowc_wc
#'
#'
#'
#' @return  an object of the same form as effect_estimate, which corresponds to the posterior mean.
get_post_effect_multfsusie <- function(G_prior, effect_estimate, list_indx_lst=NULL, lowc_wc = NULL){

  out <- list( res_uni = NULL,
               res_f  = NULL)
 if( !is.null(effect_estimate$res_uni)){
   out$res_uni$Bhat  <- do.call( cbind, lapply( 1:length(G_prior$G_prior_u), function( k) get_post_mean_u(G_prior$G_prior_u[[k]],
                                                                                effect_estimate$res_uni$Bhat[,k],
                                                                                effect_estimate$res_uni$Shat[,k]
                                                                                )
                                             )
   )
   out$res_uni$Shat  <- do.call( cbind, lapply( 1:length(G_prior$G_prior_u), function( k) get_post_sd_u(G_prior$G_prior_u[[k]],
                                                                                                          effect_estimate$res_uni$Bhat[,k],
                                                                                                          effect_estimate$res_uni$Shat[,k]
                                                                                                )
                                                )
                                    )
 }
  if( !is.null(effect_estimate$res_f)){
   out$res_f <- lapply( 1:length(G_prior$G_prior_f), function(k)list_post_mean_sd(G_prior$G_prior_f[[k]],
                                                               effect_estimate$res_f[[k]]$Bhat,
                                                               effect_estimate$res_f[[k]]$Shat,
                                                               list_indx_lst[[k]],
                                                               lowc_wc= lowc_wc)
                        )
  }
  return( out)
}





get_post_F  <- function(multfsusie.obj,l, ...)
  UseMethod("get_post_F")



#' @rdname get_post_F
#'
#' @method get_post_F multfsusie
#'
#' @export get_post_F.multfsusie
#'
#' @export
#'
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






#' @rdname get_post_F2
#'
#' @method get_post_F2 multfsusie
#'
#' @export get_post_F2.multfsusie
#'
#' @export
#'
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
                                multfsusie.obj$alpha[[l]] * (multfsusie.obj$fitted_uni2[[l]]^2+multfsusie.obj$fitted_uni[[l]]^2    )
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
      out_u <-   multfsusie.obj$alpha[[l]] * (multfsusie.obj$fitted_uni2[[l]]^2+multfsusie.obj$fitted_uni[[l]]^2    )

    }else{
      out_u <- NULL
    }
  }

  out_post <- list(post_uni_sd2 = out_u,
                   post_f_sd2    = out_f
  )
  return(out_post)
}



#' @rdname get_ER2
#'
#' @method get_ER2 multfsusie
#'
#' @export get_ER2.multfsusie
#'
#' @export

get_ER2.multfsusie = function (  multfsusie.obj,Y, X) {
  postF <- get_post_F(multfsusie.obj )# J by N matrix
  #Xr_L = t(X%*% postF)
  postF2 <- get_post_F2(multfsusie.obj ) # Posterior second moment.

  ER2 <-  list()
  if(! is.null(Y$Y_u)){

    ER2$uni <-  do.call( c,
                          lapply(1:ncol( Y$Y_u),
                                function(k) sum((Y$Y_u[,k] - X%*%postF$post_uni[,k] )^2) -sum( postF$post_uni_sd2[,k]^2) +sum( postF2$post_uni_sd2[,k])
                          )
    )


  }else
  {
    ER2$uni <- NULL
  }
  if( !is.null(Y$Y_f))
  {
    ER2$f <-  do.call( c,
                         lapply(1:length( Y$Y_f),
                                function(k) sum((Y$Y_f[[k]] - X%*%postF$post_f[[k]])^2)  -sum(postF$post_f[[k]]^2) + sum( postF2$post_f_sd2 [[k]])
                         )
    )

  }else{
    ER2$f <- NULL
  }
  return(ER2)
}



#' @title Update  multfsusie via greedy search or backfit
#'
#' @param multfsusie.obj a multfsusie object defined by \code{\link{init_multfsusie_obj}} function
#'
#' @return multfsusie object
#'
#' @export
#'
#'
greedy_backfit  <-  function(multfsusie.obj, verbose,cov_lev,X,min.purity, ...  )
  UseMethod("greedy_backfit")

#' @rdname greedy_backfit
#'
#' @method greedy_backfit multfsusie
#'
#' @export greedy_backfit.multfsusie
#'
#' @export
#'

greedy_backfit.susiF <-  function(multfsusie.obj,verbose,cov_lev,X,min.purity, ...  )
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
      print( paste( "Discarding ", length(dummy.cs), " effects"))

      multfsusie.obj <- discard_cs(multfsusie.obj,
                                   cs= dummy.cs,
                                   out_prep= FALSE
      )
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
  }

  if( length(dummy.cs)==0& !( multfsusie.obj$greedy))
  {
    multfsusie.obj$backfit <- FALSE
  }

  if(!(multfsusie.obj$greedy )&!(multfsusie.obj$backfit ) ){
    if(verbose){
      print( "Greedy search and backfitting done")
    }
    multfsusie.obj$greedy_backfit_update <- FALSE

    return(multfsusie.obj)
  }
  if(multfsusie.obj$greedy & (length(dummy.cs)==0)){

    tt <- multfsusie.obj$L_max -multfsusie.obj$L
    temp <- min( ifelse(tt>0,tt,0 ) , 7)

    if(temp==0){
      if(verbose){
        print( "Greedy search and backfitting done")
      }
      multfsusie.obj$greedy_backfit_update <- FALSE
      multfsusie.obj$backfit <- FALSE
      multfsusie.obj$greedy <- FALSE
      return(multfsusie.obj)
    }


    if(verbose){
      print( paste( "Adding ", temp, " extra effects"))
    }

    multfsusie.obj <- expand_multfsusie_obj(multfsusie.obj,L_extra = temp)
    return(multfsusie.obj)
  }

}

#' @title formatting function for output of posterior quantities
#' @description formatting function for output of posterior quantities
#' @param G_prior a mixutre_per_scale prior
#' @param  Bhat matrix of estimated mean
#' @param Shat matrix of estimated sd
#'

list_post_mean_sd <- function(G_prior, Bhat,Shat,  indx_lst, lowc_wc=NULL)
{
  out <- list (Bhat= post_mat_mean( G_prior ,
                                    Bhat,
                                    Shat,
                                    indx_lst,
                                    lowc_wc=NULL),
               Shat= post_mat_sd(   G_prior ,
                                    Bhat,
                                    Shat,
                                    indx_lst,
                                    lowc_wc=NULL)
  )
  return(out)

}




out_prep.multfsusie <- function(multfsusie.obj,Y,X,list_indx_lst,filter.cs )
{
  multfsusie.obj <-  update_cal_pip(multfsusie.obj)
  multfsusie.obj <-  update_cal_cs(multfsusie.obj)

  if(filter.cs)
  {
    multfsusie.obj <- check_cs(multfsusie.obj,min.purity=0.5,X=X)
  }
  return( multfsusie.obj)
}



pred_partial_u <- function( multfsusie.obj, l, X )
{
  tbeta <-multfsusie.obj$alpha[[l]] *multfsusie.obj$fitted_uni[[l]]
  pred_l   <- X%*% ( tbeta)
  return(pred_l)
}




#' @title Update alpha_hist   multfsusie object
#'
#' @param multfsusie.obj a multfsusie object defined by \code{\link{init_multfsusie_obj}} function
#'
#' @return multfsusie object
#'
#' @export
#'

update_alpha_hist  <-  function(multfsusie.obj, ... )
  UseMethod("update_alpha_hist")


#' @rdname update_alpha
#'
#' @method update_alpha_hist multfsusie
#'
#' @export update_alpha_hist.multfsusie
#'
#' @export
#'
update_alpha_hist.multfsusie <-  function(multfsusie.obj , ... )
{

  multfsusie.obj$alpha_hist[[ (length(multfsusie.obj$alpha_hist)+1)  ]] <- multfsusie.obj$alpha
  return( multfsusie.obj)
}


#'@title Update  multfsusie object using the output of EM_pi
#'
#' @param multfsusie.obj a multfsusie object defined by \code{\link{init_multfsusie_obj }} function
#'
#' @param l integer larger or equal to 1. Corresponds to the effect to be accessed
#'
#' @param EM_pi an object of the class "EM_pi" generated by the function \code\link{EM_pi_multfsusie}}
#' @param effect_estimate a list of marginal association generated by \code{cal_Bhat_Shat_multsusie}
#' @param list_indx_lst list of list generated by \code{\link{gen_wavelet_indx}} for the given level of resolution
#'
#' @return multfsusie object
#'
#' @export

update_multfsusie   <- function(multfsusie.obj, l, EM_pi, effect_estimate, list_indx_lst,  ...)
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
                                             lowc_wc         = lowc_wc)

  if(!is.null(post_effect$res_uni)){
    multfsusie.obj$ fitted_uni[[l]]        <-   post_effect$res_u$Bhat
    multfsusie.obj$ fitted_uni2[[l]]       <-   post_effect$res_u$Shat
  }

  if(!is.null(post_effect$res_f)){
    multfsusie.obj$fitted_wc[[l]]   <-  lapply(1: length( post_effect$res_f) ,
                                               function(k) post_effect$res_f[[k]]$Bhat)

    multfsusie.obj$fitted_wc2[[l]]  <- lapply(1: length( post_effect$res_f) ,
                                              function(k) (post_effect$res_f[[k]]$Shat^2))

  }

  new_alpha      <- susiF.alpha::cal_zeta (  EM_pi$lBF)
  multfsusie.obj <- update_alpha(multfsusie.obj, l, new_alpha)
  multfsusie.obj <- update_lBF(multfsusie.obj, l, EM_pi$lBF)

  return(multfsusie.obj)
}


#' @rdname update_pi
#'
#' @method update_pi multfsusie
#'
#' @export update_pi.multfsusie
#'
#' @export
#'
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





#' @rdname update_KL
#'
#' @method update_KL multfsusie
#'
#' @export update_KL.multfsusie
#'
#' @export
#'

update_KL.multfsusie <- function(multfsusie.obj, Y, X , list_indx_lst, ...)
{
  multfsusie.obj$KL <-  do.call(c,lapply(1:multfsusie.obj$L,
                                         FUN=function(l)
                                           cal_KL_l(multfsusie.obj, l, Y, X, list_indx_lst)))
  return( multfsusie.obj)
}



#'@title Update multfsusie log Bayes factor
#'
#'@param multfsusie.obj a multfsusie object defined by \code{\link{init_multfsusie_obj}} function
#'@param  ELBO new ELBO value
#'@return multfsusie object
#'@export

update_ELBO  <- function    (multfsusie.obj,ELBO , ...)
  UseMethod("update_ELBO")


#' @rdname update_ELBO
#'
#' @method update_ELBO multfsusie
#'
#' @export update_ELBO.multfsusie
#'
#' @export
#'

update_ELBO.multfsusie <- function    (multfsusie.obj,ELBO, ...)
{

  multfsusie.obj$ELBO <- c(multfsusie.obj$ELBO,ELBO)
  return(multfsusie.obj)
}

update_residual_variance  <- function(multfsusie.obj,sigma2, ...)
  UseMethod("update_residual_variance")



#' @rdname update_residual_variance
#'
#' @method update_residual_variance multfsusie
#'
#' @export update_residual_variance.multfsusie
#'
#' @export
#'

update_residual_variance.multfsusie <- function(multfsusie.obj,sigma2)
{
  multfsusie.obj$sigma2 <- sigma2
  return(multfsusie.obj)
}




#'@title Update multfsusie by computing PiP
#'
#'@param multfsusie.obj a multfsusie object defined by \code{\link{init_multfsusie_obj}} function
#'@return multfsusie object
#'@export

update_cal_pip  <- function (multfsusie.obj, ...)
  UseMethod("update_cal_pip")

#' @rdname update_cal_pip
#'
#' @method update_cal_pip multfsusie
#'
#' @export update_cal_pip.multfsusie
#'
#' @export
#'

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




#'@title Update multfsusie by computing credible sets
#'
#' @param multfsusie.obj a multfsusie object defined by \code{\link{init_multfsusie_obj}} function
#'
#' @param cov_lev numeric between 0 and 1, corresponding to the expected level of coverage of the cs if not specified set to 0.95
#'
#' @return multfsusie object
#'
#' @export

update_cal_cs  <- function(multfsusie.obj, cov_lev=0.95, ...)
  UseMethod("update_cal_cs")

#' @rdname update_cal_cs
#'
#' @method update_cal_cs multfsusie
#'
#' @export update_cal_cs.multfsusie
#'
#' @export
#'

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
#'
#' @param multfsusie.obj a multfsusie object defined by \code{\link{init_multfsusie_obj}} function
#'
#' @param Y  data
#' @param X matrix of size N by p
#'
#' @param list_indx_lst list generated by gen_wavelet_indx for the given level of resolution
#'
#' @param filter.cs logical, if TRUE filter the credible set (removing low purity cs and cs with estimated prior equal to 0)
#'
out_prep <- function(multfsusie.obj,Y,X,list_indx_lst,filter.cs, ...)
  UseMethod("out_prep")

#' @rdname out_prep
#'
#' @method out_prep multfsusie.obj
#'
#' @export out_prep.multfsusie.obj
#'
#' @export
#'






#' @rdname update_alpha
#'
#' @method update_alpha multfsusie
#'
#' @export update_alpha.multfsusie
#'
#' @export
#'
update_alpha.multfsusie <-  function(multfsusie.obj, l, alpha, ... )
{
  multfsusie.obj$alpha[[l]] <- alpha
  multfsusie.obj$alpha_hist[[ (length(multfsusie.obj$alpha_hist)+1)  ]] <- alpha
  return( multfsusie.obj)
}

#' @rdname update_lBF
#'
#' @method update_lBF multfsusie
#'
#' @export update_lBF.multfsusie
#'
#' @export
#'

update_lBF.multfsusie <- function    (multfsusie.obj,l, lBF, ...)
{
  if(l> multfsusie.obj$L)
  {
    stop("Error: trying to update more effects than the number of specified effect")
  }

  multfsusie.obj$lBF[[l]] <- lBF
  return(multfsusie.obj)
}





#'
#' @title Return which credible sets are  dummy
#'
#' @param multfsusie.obj a susif object defined by \code{\link{init_susiF_obj}} function
#' @param min.purity minimal purity within a CS
#' @param X matrix of covariate
#'
#' @return a list of index corresponding the the dummy effect
#'
#' @export
#'
#'
which_dummy_cs <- function(multfsusie.obj, min.purity=0.5,X,...)
  UseMethod("which_dummy_cs")

which_dummy_cs.multfsusie  <- function(multfsusie.obj, min.purity =0.5, X){
  dummy.cs<- c()

  if(is.null(multfsusie.obj$G_prior$G_prior_u)){
    for (l in 1:multfsusie.obj$L )
    {

      if (length(multfsusie.obj$cs[[l]])==1)
      {

        if(  mean( sapply(
          sapply(multfsusie.obj$est_pi[[l]]$est_pi_f,"[[",1)
          ,"[[",1)
        )==1
        ){# check if the estimated prior is exactly 0

          dummy.cs<-  c( dummy.cs,l)
        }

      }else{

        if( min(cor( X[,multfsusie.obj$cs[[l]]])) <  min.purity){#check if the purity of cs l is lower that min.purity

          dummy.cs<-  c( dummy.cs,l)

        }else{
          if(  mean( sapply(
            sapply(multfsusie.obj$est_pi[[l]]$est_pi_f,"[[",1)
            ,"[[",1)
          )==1
          ){
            dummy.cs<-  c( dummy.cs,l)
          }

        }
      }

    }
  }




  if(is.null(multfsusie.obj$G_prior$G_prior_f))
  {
    for (l in 1:multfsusie.obj$L )
    {

      if (length(multfsusie.obj$cs[[l]])==1)
      {

        if(  mean( sapply(
          sapply(multfsusie.obj$est_pi[[l]]$est_pi_u,"[[",1)
          ,"[[",1)
        )==1
        ){# check if the estimated prior is exactly 0

          dummy.cs<-  c( dummy.cs,l)
        }

      }else{

        if( min(cor( X[,multfsusie.obj$cs[[l]]])) <  min.purity){#check if the purity of cs l is lower that min.purity

          dummy.cs<-  c( dummy.cs,l)

        }else{
          if( mean( sapply(
            sapply(multfsusie.obj$est_pi[[l]]$est_pi_u,"[[",1)
            ,"[[",1)
          )==1){
            dummy.cs<-  c( dummy.cs,l)
          }

        }
      }

    }

  }


  if ((!is.null(multfsusie.obj$G_prior$G_prior_f))& (!is.null(multfsusie.obj$G_prior$G_prior_u))){
    for (l in 1:multfsusie.obj$L )
    {

      if (length(multfsusie.obj$cs[[l]])==1)
      {

        m_u <-mean( sapply(
          sapply(multfsusie.obj$est_pi[[l]]$est_pi_u,"[[",1)
          ,"[[",1)
        )
        m_f <- mean( sapply(
          sapply(multfsusie.obj$est_pi[[l]]$est_pi_u,"[[",1)
          ,"[[",1)
        )

        if( mean(c(m_u,m_f))==1){# check if the estimated prior is exactly 0

          dummy.cs<-  c( dummy.cs,l)
        }

      }else{

        if( min(cor( X[,multfsusie.obj$cs[[l]]])) <  min.purity){#check if the purity of cs l is lower that min.purity

          dummy.cs<-  c( dummy.cs,l)

        }else{
          m_u <-mean( sapply(
            sapply(multfsusie.obj$est_pi[[l]]$est_pi_u,"[[",1)
            ,"[[",1)
          )
          m_f <- mean( sapply(
            sapply(multfsusie.obj$est_pi[[l]]$est_pi_u,"[[",1)
            ,"[[",1)
          )

          if( mean(c(m_u,m_f))==1){
            dummy.cs<-  c( dummy.cs,l)
          }

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
