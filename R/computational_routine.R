

#' @title Compute Bhat / Shat across all multfSuSiE modalities
#'
#' @param Y          named list with optional Y_u (N x K_u, univariate traits)
#'                   and optional Y_f (length-M list of N x T_m wavelet
#'                   coefficient matrices).
#' @param X          N x p predictor matrix.
#' @param sigma2     multfsusie.obj$sigma2 — a list with components
#'                     $sd_u  (length K_u, per univariate trait), and/or
#'                     $sd_f  (length M scalars per modality, or a length-M
#'                            list of length-T_m vectors if you have moved
#'                            to per-position residual variance).
#' @param low_trait  list of low-count masks: $low_u (vec of trait indices to
#'                   mask) and $low_wc (length-M list of position-index vecs).
#' @param ind_analysis  named list with $idx_u (length K_u of row vectors)
#'                      and $idx_f (length M of row vectors); missing or NULL
#'                      means "use all rows".
#' @param v1, list_indx_lst   ignored — kept for backward signature
#'   compatibility.
#'
#' @return list(res_u = list(Bhat, Shat) | NULL,
#'              res_f = list of M lists(Bhat, Shat) | NULL)
#' @export
## ============================================================================
## cal_Bhat_Shat_multfsusie  —  zero-pad-and-correct version
##
## Designed for M up to many hundreds of modalities, each with its own row
## mask.  Replaces M small crossprod calls with ONE large one, regardless of
## how heterogeneous the missing-data pattern is across modalities.
##
## Algorithm (functional block, full data — easy case):
##   d <- attr(X, "d")                           # cached column-norms-squared
##   Y_stack <- do.call(cbind, Y_f)              # N x sum(T_m)
##   B_stack <- crossprod(X, Y_stack) / d        # ONE BLAS call
##   slice B_stack into per-modality (Bhat_m, Shat_m)
##
## Algorithm (functional block, with per-modality row masks):
##   pad each Y_f[[m]] with zeros at the rows NOT in idx_f[[m]]:
##       Y_pad_m[i, ] = Y_f[[m]][i, ] if i in idx_f[[m]] else 0
##   Y_stack <- do.call(cbind, Y_pad)            # still N x sum(T_m)
##   B_stack <- crossprod(X, Y_stack)            # ONE BLAS call, no /d yet
##
##   for each modality m:
##     d_m <- attr(X,"d") - colSums(X[!idx_m, , drop=FALSE]^2)   # incremental
##     Bhat_m <- B_stack[, cols_m] / d_m                          # per-m
##     Shat_m <- (1/sqrt(d_m)) %o% sqrt(sigma2_m)                 # tcrossprod
##
## The d_m correction is O(|missing| * p) per modality.  When missingness is
## sparse (typical), this is negligible.  When it is dense, fall back to a
## direct colSums(X[idx_m, ]^2) per modality (also O(N*p) per m, but you'd
## have paid that anyway in the per-m approach — only the BLAS overhead is
## now amortised).
##
## Memory note:
##   Y_stack is an extra N x sum(T_m) double array.  For N=1000, T_m=128,
##   M=500 that's ~4 GB.  If memory is constrained, use chunk_size below to
##   process modalities in chunks of K -> M/K BLAS calls instead of M.
##   K=50 typically keeps the working set under 500 MB and still gives a
##   ~10x speedup over per-m.
## ============================================================================

cal_Bhat_Shat_multfsusie <- function(Y, X, sigma2,
                                     low_trait     = NULL,
                                     ind_analysis  = NULL,
                                     chunk_size    = Inf,    # set to e.g. 50 if memory tight
                                     v1            = NULL,   # ignored
                                     list_indx_lst = NULL,   # ignored
                                     ...) {

  if (is.null(sigma2))
    stop("cal_Bhat_Shat_multfsusie: pass multfsusie.obj$sigma2 as 'sigma2'.")
  if (!is.list(sigma2))
    stop("cal_Bhat_Shat_multfsusie: 'sigma2' must be a list with $sd_u/$sd_f. ",
         "Did you pass legacy `v1` positionally?")

  has_ind   <- !is.null(ind_analysis)
  has_idx_u <- has_ind && !is.null(ind_analysis$idx_u)
  has_idx_f <- has_ind && !is.null(ind_analysis$idx_f)

  ## --- Univariate block ----------------------------------------------------
  if (is.null(Y$Y_u)) {
    res_u <- NULL
  } else {
    if (is.null(sigma2$sd_u))
      stop("cal_Bhat_Shat_multfsusie: Y$Y_u present but sigma2$sd_u is NULL.")
    res_u <- fsusieR:::cal_Bhat_Shat(
      Y            = Y$Y_u,
      X            = X,
      sigma2       = sigma2$sd_u,
      lowc_wc      = low_trait$low_u,
      ind_analysis = if (has_idx_u) ind_analysis$idx_u else NULL
    )
  }

  ## --- Functional block ----------------------------------------------------
  if (is.null(Y$Y_f)) {
    res_f <- NULL
  } else {
    if (is.null(sigma2$sd_f))
      stop("cal_Bhat_Shat_multfsusie: Y$Y_f present but sigma2$sd_f is NULL.")

    res_f <- .stacked_functional_block(
      X         = X,
      Y_f       = Y$Y_f,
      sigma2_f  = sigma2$sd_f,
      idx_f     = if (has_idx_f) ind_analysis$idx_f else NULL,
      low_wc    = low_trait$low_wc,
      chunk_size = chunk_size
    )
  }

  list(res_u = res_u, res_f = res_f)
}


## ----------------------------------------------------------------------------
## Stacked driver — one BLAS call per chunk of modalities.
## ----------------------------------------------------------------------------
.stacked_functional_block <- function(X, Y_f, sigma2_f, idx_f, low_wc,
                                      chunk_size = Inf) {

  M  <- length(Y_f)
  N  <- nrow(X)
  P  <- ncol(X)
  Tm <- vapply(Y_f, ncol, integer(1))           # length M

  ## Cached full-data column norms (set by fsusieR::colScale)
  d_full <- attr(X, "d")
  if (is.null(d_full)) d_full <- .colSums(X * X, N, P)

  sd_f_is_list <- is.list(sigma2_f)
  get_sigma_m  <- function(m) if (sd_f_is_list) sigma2_f[[m]] else sigma2_f[m]
  get_idx_m    <- function(m) if (is.null(idx_f)) NULL else idx_f[[m]]

  res_f <- vector("list", M)

  ## Chunk modalities to bound peak memory of Y_stack.
  chunk_size  <- max(1L, min(M, as.integer(chunk_size)))
  chunk_starts <- seq.int(1L, M, by = chunk_size)

  for (cs in chunk_starts) {
    ce  <- min(cs + chunk_size - 1L, M)
    grp <- cs:ce

    ## ---- Build Y_stack for this chunk: zero-pad missing rows ------------
    ##
    ## We allocate Y_stack as N x sum(Tm[grp]) and copy each Y_f[[m]] into
    ## its column slice, zeroing rows that are missing for that modality.
    Tm_grp     <- Tm[grp]
    col_end    <- cumsum(Tm_grp)
    col_start  <- c(1L, head(col_end + 1L, -1L))
    T_total    <- sum(Tm_grp)
    Y_stack    <- matrix(0, nrow = N, ncol = T_total)

    for (i in seq_along(grp)) {
      m  <- grp[i]
      cs_cols <- col_start[i]:col_end[i]
      idx_m   <- get_idx_m(m)
      if (is.null(idx_m)) {
        ## All rows used.  Direct copy.  If Y_f[[m]] still has NAs, scrub.
        Yi <- Y_f[[m]]
        if (anyNA(Yi)) Yi[is.na(Yi)] <- 0
        Y_stack[, cs_cols] <- Yi
      } else {
        ## Only observed rows contribute; the rest stay zero.
        Y_stack[idx_m, cs_cols] <- Y_f[[m]][idx_m, , drop = FALSE]
        ## (If Y_f[[m]] has NAs at observed rows, that's a data bug —
        ## let it propagate so the caller notices.)
      }
    }

    ## ---- ONE BLAS call for the whole chunk ------------------------------
    ## We hold off on dividing by d here because d is per-modality.
    B_stack <- crossprod(X, Y_stack)                      # P x T_total

    ## ---- Slice + per-modality d_m + finalize ----------------------------
    for (i in seq_along(grp)) {
      m       <- grp[i]
      cs_cols <- col_start[i]:col_end[i]
      idx_m   <- get_idx_m(m)

      ## Per-modality column norms.
      ## Sparse-missingness shortcut: d_m = d_full - sum over missing rows.
      if (is.null(idx_m)) {
        d_m <- d_full
      } else {
        n_obs   <- length(idx_m)
        n_miss  <- N - n_obs
        if (n_miss == 0L) {
          d_m <- d_full
        } else if (n_miss * 4L < N) {
          ## Fewer missing than 1/4 of rows — incremental subtraction is
          ## cheaper than rebuilding from idx_m.
          miss_idx <- setdiff(seq_len(N), idx_m)
          Xm <- X[miss_idx, , drop = FALSE]
          d_m <- d_full - .colSums(Xm * Xm, n_miss, P)
        } else {
          ## Lots of missing — direct recomputation on observed rows.
          Xo <- X[idx_m, , drop = FALSE]
          d_m <- .colSums(Xo * Xo, n_obs, P)
        }
      }

      ## Numerical guard: any d_m[k] = 0 means column k is zero on the
      ## modality's observed rows.  Mark via Inf Shat so BF contribution -> 0.
      bad <- d_m <= 0
      if (any(bad)) d_m[bad] <- 1   # avoid /0; we'll overwrite Shat below

      Bhat_m <- B_stack[, cs_cols, drop = FALSE] / d_m

      sigma2_m <- rep_len(as.numeric(get_sigma_m(m)), Tm_grp[i])
      Shat_m   <- tcrossprod(1 / sqrt(d_m), sqrt(sigma2_m))   # P x T_m

      if (any(bad)) {
        Bhat_m[bad, ] <- 0
        Shat_m[bad, ] <- 1
      }

      lw <- if (is.null(low_wc)) NULL else low_wc[[m]]
      if (!is.null(lw)) {
        Bhat_m[, lw] <- 0
        Shat_m[, lw] <- 1
      }

      Shat_m[Shat_m < 1e-32] <- 1e-32

      res_f[[m]] <- list(Bhat = Bhat_m, Shat = Shat_m)
    }
  }

  res_f
}





#' @title Compute conditional local false sign rate
#
#' @description Compute conditional local false sign rate
#'
#' @param G_prior multfsusie_prior
#
#' @param effect_estimate output of cal_Bhat_Shat_multfsusie
#'
#' @param list_indx_lst list of list generated by \code{\link{gen_wavelet_indx}} for the given level of resolution
#'
#' @param \ldots other arguments
#'
#' @return esitmated conditional lfsr
#
#' @export


cal_clfsr <- function (G_prior, effect_estimate, list_indx_lst,...)
  UseMethod("cal_clfsr")


#' @rdname cal_clfsr
#
#' @method cal_clfsr multfsusie_prior
#
#' @export cal_clfsr.multfsusie_prior
#
#
#
#' @importFrom ashr set_data
#' @importFrom ashr get_fitted_g
#' @importFrom fsusieR cal_clfsr.mixture_normal_per_scale
#  @importFrom ashr calc_lfsr
#' @export
#

cal_clfsr.multfsusie_prior <- function(G_prior ,
                                       effect_estimate,
                                       list_indx_lst,...){

  if( ! is.null(effect_estimate$res_f)){
    clfsr_wc <- NULL


    #TODO the line below take to much memory
    # this need to be cleaned
    #clfsr_wc <- lapply(1: length(effect_estimate$res_f),
    #                    function(k){
    #
    #                     fsusieR::cal_clfsr (
    #                       G_prior  = G_prior$G_prior_f[[k]],
    #                       Bhat     = effect_estimate$res_f[[k]]$Bhat,
    #                       Shat     = effect_estimate$res_f[[k]]$Shat,
    #                       indx_lst = list_indx_lst[[k]]
    #                     )
    #
    #                  }
    #)
  }else{
    clfsr_wc <- NULL
  }

  if( ! is.null(effect_estimate$res_u)){
    clfsr_u <- do.call(rbind,
                         lapply(1:ncol(effect_estimate$res_u$Bhat),
                                function(k){
                                  m <- G_prior$G_prior_u[[k]] [[1]]

                                  data_ash <-  ashr::set_data(
                                    effect_estimate$res_u$Bhat[,k],
                                    effect_estimate$res_u$Shat[,k])

                                  ashr:::calc_lfsr( m ,data_ash)
                                }
                         )
    )

  }else{
    clfsr_u <- NULL
  }

  clfsr_mult <- list( clfsr_wc   = clfsr_wc,
                      clfsr_u = clfsr_u)
  out <-  clfsr_mult
  return(out)

}











# @title Compute Log-Bayes Factor for univariate regression with ash prior
#
# @description Compute Log-Bayes Factor
#
# @param G_prior ash object
#
# @param Bhat p numerical vector of regression coefficients;
#
# @param Shat p numerical of standard errors;
# @return  The log-Bayes factor for each covariate.
#
# @export

log_BFu <- function (G_prior, Bhat, Shat,low_u=FALSE,df=NULL, ...) {


  Shat[ Shat<=0 ] <- 1e-32

  if( is.null(df)){
    if( low_u){
      out <- rep(0, length(Bhat))
    }else{
      tt   <- rep(0,length(Shat))
      pi_k <- G_prior$fitted_g$pi
      sd_k <- G_prior$fitted_g$sd
      # Speed Gain: could potential skip the one that are exactly zero.
      for (o in 1:length(G_prior$fitted_g$pi)){
        tt <- tt + pi_k[o] * stats::dnorm(Bhat ,sd = sqrt(sd_k[o]^2 + Shat ^2))
      }

      out <-  (log(tt) - stats::dnorm(Bhat ,sd = Shat ,log = TRUE))

    }
  }else{
    if( low_u){
      out <- rep(0, length(Bhat))
    }else{
      tt   <- rep(0,length(Shat))
      pi_k <- G_prior$fitted_g$pi
      sd_k <- G_prior$fitted_g$sd
      # Speed Gain: could potential skip the one that are exactly zero.
      for (o in 1:length(G_prior$fitted_g$pi)){
        tt <- tt + pi_k[o] *LaplacesDemon::dstp(Bhat,tau = 1/(sd_k[o]^2 + Shat ^2), nu=df)
      }
      out <- sum(log(tt) - LaplacesDemon::dstp(Bhat ,tau = 1/Shat ^2,nu=df,log = TRUE))
    }


  }



  return(out)
}

log_BF <- function( G_prior,effect_estimate ,list_indx_lst,low_trait , df=NULL )
  UseMethod("log_BF")

# @title Compute Log-Bayes Factor for a multiple f susie regression model
# @description Compute Log-Bayes Factor
#
# @param G_prior a multfsusie_prior
#
# @param effect_estimate regression coefficients generated by \link{\code{cal_Bhat_Shat_multfsusie}}
#
# @param  list_indx_lst List of lists generated by \code{\link{gen_wavelet_indx}}
#   for the given level of resolution
# @return  The log-Bayes factor for each covariate.
#
#' @export
#' @keywords internal
log_BF.multfsusie_prior <- function( G_prior,
                                     effect_estimate ,
                                     list_indx_lst,
                                     low_trait,
                                     df=NULL)
{
  if(is.null(df)){
    if( is.null(G_prior$G_prior_u)){
      u_logBF <-  matrix(rep(0,nrow(effect_estimate$res_f[[1]]$Bhat  )), nrow=1)
    }else{
      u_logBF <-  lapply(1:ncol(effect_estimate$res_u$Bhat),
                         function(k)
                           log_BFu(G_prior = G_prior$G_prior_u[[k]],
                                   Bhat    =  effect_estimate$res_u$Bhat[,k] ,
                                   Shat    =  effect_estimate$res_u$Shat[,k],
                                   low_u   =  ifelse(k%in%low_trait$low_u, TRUE,FALSE)
                           )
      )
      u_logBF <-  do.call(rbind, u_logBF)
    }
    if(is.null(G_prior$G_prior_f)){
      f_logBF <-  matrix(rep(0,nrow(effect_estimate$res_u[[1]] )), nrow=1)
    }else{
      f_logBF <- lapply( 1: length(G_prior$G_prior_f) ,function( k)
        fsusieR::log_BF(G_prior  = G_prior$G_prior_f[[k]],
                            Bhat     = effect_estimate$res_f[[k]]$Bhat,
                            Shat     = effect_estimate$res_f[[k]]$Shat,
                            indx_lst = list_indx_lst[[k]],
                            lowc_wc  = low_trait$lowc_wc[[k]]
        )
      )
      f_logBF <- do.call(rbind, f_logBF)
    }

  }else{
    ### Work here ------
    if( is.null(G_prior$G_prior_u)){
      u_logBF <- matrix(rep(0,nrow(effect_estimate$res_f[[1]]$Bhat  )), nrow=1)
    }else{
     # print( "moderated BF u ")
      u_logBF <-  lapply(1:ncol(effect_estimate$res_u$Bhat),
                         function(k)
                           log_BFu(G_prior =  G_prior$G_prior_u[[k]],
                                   Bhat    =  effect_estimate$res_u$Bhat[,k] ,
                                   Shat    =  effect_estimate$res_u$Shat[,k],
                                   low_u   =  ifelse(k%in%low_trait$low_u, TRUE,FALSE),
                                   df      =  df$Y_u[k]
                           )
      )
      u_logBF <-  do.call(rbind, u_logBF)
    }
    if(is.null(G_prior$G_prior_f)){
      f_logBF <- matrix(rep(0,nrow(effect_estimate$res_u[[1]] )), nrow=1)
    }else{
      #print( "moderated BF f ")
      f_logBF <- lapply( 1: length(G_prior$G_prior_f) ,function( k)
        fsusieR::log_BF(G_prior  = G_prior$G_prior_f[[k]],
                            Bhat     = effect_estimate$res_f[[k]]$Bhat,
                            Shat     = effect_estimate$res_f[[k]]$Shat,
                            indx_lst = list_indx_lst[[k]],
                            lowc_wc  = low_trait$lowc_wc[[k]],
                            df       = df$Y_f[k]
        )
      )
      f_logBF <-do.call(rbind, f_logBF)
    }



  }

  out <- list(  f_logBF=f_logBF,
                u_logBF=u_logBF)
  return(out)

}







# @title Compute posterior mean for univariate regression
# @description Compute posterior mean for univariate regression
# @param Bhat  a vector of mean estimate
# @param Bhat  a vector of sd estimate
# @param low_u logical indicate if the trait as critically low spread
get_post_mean_u <- function(G_prior, Bhat, Shat, low_u=FALSE)
{
  if(low_u){
    return(rep( 0, length(Bhat)))
  }else{
    data <-  ashr::set_data(Bhat  ,Shat  )
    return(ashr::postmean(ashr::get_fitted_g(G_prior),data))
  }

}


# @title Compute posterior sd for univariate regression
# @description Compute posterior sd for univariate regression
# @param Bhat  a vector of mean estimate
# @param Bhat  a vector of sd estimate
# @param low_u logical indicate if the trait as critically low spread
get_post_sd_u <- function(G_prior, Bhat, Shat, low_u=FALSE)
{
  if(low_u){
    return(rep( 1, length(Bhat)))
  }else{
  data <-  ashr::set_data(Bhat  ,Shat  )
  return(ashr::postsd(ashr::get_fitted_g(G_prior),data))
  }
}


#' @title  Compute Residual variance
#' @description  see title
#' @param multfsusie.obj a multfsusie object
#' @param Y observed response data
#' @param X observed covariates
#' @export
#' @keywords internal
estimate_residual_variance <- function(multfsusie.obj,Y,X,... )
  UseMethod("estimate_residual_variance")



#' @rdname estimate_residual_variance
#
#' @method estimate_residual_variance multfsusie
#
#' @export estimate_residual_variance.multfsusie
#
#' @export
#' @keywords internal
estimate_residual_variance.multfsusie <- function(multfsusie.obj,Y,X, ind_analysis, ... )
{
  if (missing(ind_analysis)){
    R2 <- get_ER2( multfsusie.obj, Y, X)
    est_sd2 <-  list()
    if(!is.null(R2$uni))
    {
      est_sd2$sd_u <-  R2$uni/nrow(Y$Y_u)
    }
    if(!is.null(R2$f)){
      n <- rep(nrow(Y$Y_f[[1]]), length(Y$Y_f) )
      t <- do.call( c, lapply(1: length(Y$Y_f), function(k) ncol(Y$Y_f[[k]] ) ))
      est_sd2$sd_f <- R2$f / (n*t)
    }

  }else{
    R2 <- get_ER2( multfsusie.obj, Y, X, ind_analysis)
    est_sd2 <-  list()
    if(!is.null(R2$uni))
    {
      est_sd2$sd_u <-  R2$uni/(nrow(Y$Y_u)- (  nrow(Y$Y_u) - lengths(ind_analysis$idx_u) ))# accounting for missing data points
    }
    if(!is.null(R2$f)){# accounting for missing data points
      n <-  do.call(c,  lapply( 1: length(Y$Y_f), function( k) (nrow(Y$Y_f[[k]])- ( nrow(Y$Y_f[[k]]) - length(ind_analysis$idx_f[[k]]) ))))
      t <- do.call( c, lapply(1: length(Y$Y_f), function(k) ncol(Y$Y_f[[k]] ) ))
      est_sd2$sd_f <- R2$f / (n*t)
    }

  }





  out <-  est_sd2
  return(out)
}






#@title Function to fit entry-wise lm on tensor regression
#
#@param l  wavelet coefficient index
#@param j  covariate index
#@param xi condition index
#@param Y observed tensor
#@param X observed covariate
# @export
parse_lm_fit <- function(j,l,xi, v1,Y,X)
{

  out <- fast_lm(cbind(v1,X[,j]),Y[,l,xi])
  return(c(out$be[2,1],
           sqrt(
             Rfast::cova(out$residuals)/sum(
               (X[,j]-mean(X[,j]))^2)
           )
  )
  )
}

