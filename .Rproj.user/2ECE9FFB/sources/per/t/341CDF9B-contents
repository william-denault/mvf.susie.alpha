
#' @title Initialize a multsusie object
#'
#' @param L number of non zero coefficients An L-vector containing the indices of the
#'   nonzero coefficients.
#'
#' @param G_prior prior object defined by init_prior_multsusie function
#'
#' @param Y list of matrices of outcomes
#'
#' @param X Matrix of covariates
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
init_mvfsusie_obj <- function(L, G_prior, Y,X )
{


  fitted_wc       <-  list()
  fitted_wc2      <-  list()
  alpha           <-  list()
  alpha_hist      <-  list()
  ind_fitted_func <-  list
  cs              <-  list()
  pip             <-  rep(0, dim(X)[2])
  est_pi          <-  list()
  est_sd          <-  list()
  lfsr_wc         <-  list()
  lfsr_func       <-  list()
  L               <-  L
  G_prior         <-  G_prior
  N               <-  nrow(X)[1]
  n_wac           <-  lapply(lapply(Y,dim) ,`[[`, 2)
  n_cond          <-  length(Y)
  P               <-  ncol(X)[2]
  sigma2          <-  rep(1, length(Y))
  lBF             <-  list()
  KL              <-  rep(NA,L)
  ELBO            <-  c()
  for ( l in 1:L )
  {
    fitted_wc[[l]]        <-  array(0, dim= c(dim(X)[2], ncol=dim(Y)[2], dim(Y)[3])  )
    fitted_wc2[[l]]       <-  array(0, dim= c(dim(X)[2], ncol=dim(Y)[2], dim(Y)[3])  )
    alpha [[l]]           <-  rep(0, dim(X)[2])
    cs[[l]]               <-  list()
    est_pi [[l]]          <-  get_pi_G_prior(G_prior)
    lBF[[l]]              <-  rep(NA, ncol(X))
    lfsr_wc[[l]]          <-  rep(1, ncol(Y))
    lfsr_func[[l]]        <-  rep(1, ncol(Y))
  }
  obj <- list( fitted_wc       = fitted_wc,
               fitted_wc2      = fitted_wc2,
               lfsr_wc         = lfsr_wc,
               lBF             = lBF,
               KL              = KL,
               ELBO            = ELBO,
               ind_fitted_func = ind_fitted_func,
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
               L               = L)

  class(obj) <- "mvfsusie"
  return(obj)
}
