
# @title Method to fit mash
#
# @description Method to fit mash
#
# @param Bhat  matrix of the regression coefficient (MLE)
# @param Shat  matrix of standard error
# @param data.driven logical, if TRUE use data driven covariance matrix procedure for fitting mash otherwise uses  cov_canonical method mashr. Set as TRUE by default
#
# @return a mash object
#
# @importFrom mashr mash_set_data
# @importFrom mashr mash_1by1
# @importFrom mashr get_significant_results
# @importFrom mashr cov_pca
# @importFrom mashr cov_ed
# @importFrom mashr cov_canonical
# @importFrom mashr mash
#
# @export
#


basic_mash_fit <- function (Bhat, Shat, data.driven=TRUE, verbose=FALSE)
{
  data   = mash_set_data( Bhat,  Shat)

  if( data.driven)
  {
    m.1by1 = mash_1by1(data)
    strong = get_significant_results(m.1by1,0.05)
    if( !(length(strong)==0))
    {
      U.pca  = cov_pca(data,
                       npc=length(strong),
                       subset=strong
      )
      #Problem with TPCa being NAN
      U.ed   = cov_ed(data, U.pca, subset=strong)
      U.c    = cov_canonical(data)
      m      = mash(data, c(U.c,U.ed),verbose = verbose)
    }else{
      U.c    = cov_canonical(data)
      m      = mash(data, c(U.c),verbose = verbose)
    }
  }else{
    U.c    = cov_canonical(data)
    m      = mash(data, c(U.c),verbose = verbose)
  }


  return(m)
}




# @title Regress column l and condition of Y on column j of X
#
# @description Add description here.
#
# @param Y  tensor phenotype, matrix of size N by size J by xi. The underlying algorithm uses wavelets that assume that J is of the form J^2. If J is not a power of 2, susiF internally remaps the data into a grid of length 2^J
#
# @param X matrix of size n by p in
#
# @param v1 vector of 1 of length n
#
# @return list of two tensor of size pxTx xi of 2 containing the regression coefficient and standard error
#
# @importFrom stats var
#
# @export
#


cal_Bhat_Shat_tensor  <- function(Y, X, v1)
{
  Bhat  <- list()
  Shat  <- list()
  h <- 1 # looping index
  #need to parse by variable  (most inside loop) then by wave coef then by condition
  for (xi in 1:dim(Y)[3])
  {
    for ( l in 1:dim(Y)[2])
    {
      out <-  do.call( cbind,lapply( 1:dim(X)[2], function(j) parse_lm_fit( j=j,l=l,xi=xi,v1=v1, Y=Y, X=X ) ) )
      Bhat[[h]]  <- out[1,]
      Shat[[h]]  <- out[2,]
      h <- h+1
    }
  }

  tens_Bhat <- array( do.call(c, Bhat), dim = c( ncol(X),dim(Y)[2], dim(Y)[3]) )
  tens_Shat <- array( do.call(c, Shat), dim = c( ncol(X),dim(Y)[2], dim(Y)[3]) )

  out <- list( tens_Bhat = tens_Bhat,
               tens_Shat = tens_Shat)
  return(out)
}




# @title Regress different marks of Y   on X nxp
#
# @description regression coefficients (and sd) of the column wise regression
#
# @param Y a list of two, Y_u containning a N by J data frame of univariate phneotype and Y_f a k list contains a list of functionnal phenoytpes
#
# @param X matrix of size N by P in
#
# @return a nested list list of two
#
# \item{res_u}{ list of two Bhat: matrix pxJ regression coefficient, Bhat[j,t] corresponds to regression coefficient of the t univariate phneotype
#  on X[,j]; Shat is the matrix of the corresponding standard error }
#
# \item{res_f}{ a list of k in which each element contains Bhat and Shat matrix (see description in item res_u}
#
# @export


cal_Bhat_Shat_multfsusie <- function( Y,X,v1,
                                      list_indx_lst=NULL,
                                      low_trait=NULL,
                                      ind_analysis, parallel =FALSE  )
{

  if(is.null(Y$Y_u)){
    res_u <- NULL
  }else{
    if(missing(ind_analysis)){
      res_u   <- fsusieR:::cal_Bhat_Shat(Y=Y$Y_u,
                                               X=X,
                                               v1=v1,
                                               lowc_wc=low_trait$low_u)
    #  res_u$Shat <- res_u$Shat%*%diag(sqrt(multfsusie.obj$sigma2$sd_u))
    }else{
      res_u   <- fsusieR:::cal_Bhat_Shat(Y=Y$Y_u,
                                               X=X,
                                               v1=v1,
                                               lowc_wc=low_trait$low_u,
                                               ind_analysis=ind_analysis$idx_u)
     # res_u$Shat <- res_u$Shat%*%diag(sqrt(multfsusie.obj$sigma2$sd_u))
    }

  }

  if(is.null(Y$Y_f)){
    res_f <- NULL
  }else{
    if(missing(ind_analysis)){



      if(parallel){
        res_f <- parallel::mclapply(1:length(Y$Y_f),
                                    function(k) fsusieR:::cal_Bhat_Shat(Y$Y_f[[k]],
                                                                X       = X,
                                                                v1      = v1,
                                                                resid_var = multfsusie.obj$sigma2$sd_f[k],
                                                                lowc_wc = low_trait$low_wc[[k]]),
                                   mc.cores=numCores,
                                   mc.preschedule=FALSE
                                   )
      }else{
        res_f <- lapply(1:length(Y$Y_f),
                        function(k) fsusieR:::cal_Bhat_Shat(Y$Y_f[[k]],
                                                                X       = X,
                                                                v1      = v1,
                                                                resid_var = multfsusie.obj$sigma2$sd_f[k],
                                                                lowc_wc = low_trait$low_wc[[k]])


                      )
      }
    }else{


      if(parallel){
        res_f <-  parallel::mclapply( 1:length(Y$Y_f),
                        function(k) fsusieR:::cal_Bhat_Shat(Y$Y_f[[k]],
                                                                X       = X,
                                                                v1      = v1,
                                                                lowc_wc = low_trait$low_wc[[k]],
                                                                resid_var = multfsusie.obj$sigma2$sd_f[k],
                                                                ind_analysis=ind_analysis$idx_f[[k]]),
                                    mc.cores=numCores,
                                    mc.preschedule=FALSE
                                    )
      }else{
        res_f <- lapply(1:length(Y$Y_f),
                        function(k) fsusieR:::cal_Bhat_Shat(Y$Y_f[[k]],
                                                                X       = X,
                                                                v1      = v1,
                                                                lowc_wc = low_trait$low_wc[[k]],
                                                                resid_var = multfsusie.obj$sigma2$sd_f[k],
                                                                ind_analysis=ind_analysis$idx_f[[k]])
        )
      }

    }

  }

 res  <- list( res_u = res_u,
               res_f   = res_f)

}





# @title Compute Log-Bayes Factor of a given   covariate for multivariate wavelet regression
#
# @description   Compute Log-Bayes Factor of a given   covariate for multivariate wavelet regression
#
# @param G_prior a scale specific mash prior
#
# @param tens_marg a list of tensor of marginal association generated by \code{cal_Bhat_Shat_tensor}
#
# @param s scale of interest
#
# @param j  covariate of interest
#
# @param indx_list List generated by \code{\link{gen_wavelet_indx}}
#   for the given level of resolution
#
# @return The log-Bayes factor for the considered covariate.
#
# @export
#

cal_lbf_cov_mvsusie <- function(G_prior,tens_marg,s,j,indx_lst)
{
  #sum over the level specific Bayes Factor
  out <- do.call(sum,
                 lapply ( 1:length(indx_lst),
                          function(s) cal_lbf_mvfsusie_level (G_prior,
                                                              tens_marg,
                                                              s,
                                                              j,
                                                              indx_lst
                          )
                 )
  )
  return(out)
}


# @title Compute Log-Bayes Factor of given scale and covariate for multivariate wavelet regression
#
# @description   Compute Log-Bayes Factor of given scale and covariate for multivariate wavelet regression
#
# @param G_prior a scale specific mash prior
#
# @param tens_marg a list of tensor of marginal association generated by \code{cal_Bhat_Shat_tensor}
#
# @param s scale of interest
#
# @param j  covariate of interest
#
# @param indx_list List generated by \code{\link{gen_wavelet_indx}}
#   for the given level of resolution
#
# @return The log-Bayes factor of the given scale and covariate
#
# @export
#

cal_lbf_mvfsusie_level <-  function(G_prior, tens_marg,s ,j , indx_lst)
{


  Bhat   <-   matrix(tens_marg$tens_Bhat[j,indx_lst[[s]],],ncol = dim(tens_marg$tens_Bhat)[3])
  Shat   <-   matrix(tens_marg$tens_Shat[j,indx_lst[[s]],],ncol = dim(tens_marg$tens_Bhat)[3])
  data   <-   mash_set_data( Bhat,  Shat)
  m      <-   G_prior[[s]]

  loglik      <-  sum (mash_compute_vloglik(m,data))
  null_loglik <-  mvfsusie_compute_null_loglik(Bhat,Shat)
  lbf         <-   loglik -null_loglik

  return(lbf)
}






#' @title Compute conditional local false sign rate
#
#' @description Compute conditional local false sign rate
#' @param G_prior multfsusie_prior
#
#' @param effect_estimate output of cal_Bhat_Shat_multfsusie
#' @param list_indx_lst list of list generated by \code{\link{gen_wavelet_indx}} for the given level of resolution
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
                                       list_indx_lst){

  if( ! is.null(effect_estimate$res_f)){
    clfsr_wc <- lapply(1: length(effect_estimate$res_f),
                       function(k){

                         fsusieR:::cal_clfsr (
                           G_prior  = G_prior$G_prior_f[[k]],
                           Bhat     = effect_estimate$res_f[[k]]$Bhat,
                           Shat     = effect_estimate$res_f[[k]]$Shat,
                           indx_lst = list_indx_lst[[k]]
                         )

                       }
    )
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











# @title Fit mash for a given set of wavelet coefficients
#
# @description Fit mash for a given set of wavelet coefficients
#
# @param tens_marg list of tensor of marginal association generated by \code{cal_Bhat_Shat_tensor}
#
# @param s scale of interest
#
# @param data.driven logical, if TRUE use data driven covariance matrix procedure for fitting mash otherwise uses  cov_canonical method mashr. Set as TRUE by default
#
# @param indx_list list generated by \code{\link{gen_wavelet_indx}} for the given level of resolution
#
# @return a mash object
#
# @export
#


fit_mash_level <- function(tens_marg, s, indx_lst, data.driven=TRUE,verbose=FALSE)
{



  Bhat <- cbind_3Darray( tens_marg$tens_Bhat[,indx_lst[[s]],] )
  Shat <- cbind_3Darray( tens_marg$tens_Shat[,indx_lst[[s]],] )


  m <- basic_mash_fit(Bhat, Shat, data.driven = data.driven,verbose=verbose)

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
        tt <- tt + pi_k[o] * dnorm(Bhat ,sd = sqrt(sd_k[o]^2 + Shat ^2))
      }

      out <-  (log(tt) - dnorm(Bhat ,sd = Shat ,log = TRUE))

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

    }
    out <- sum(log(tt) - LaplacesDemon::dstp(Bhat ,tau = 1/Shat ^2,nu=df,log = TRUE))

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
# @export
log_BF.multfsusie_prior <- function( G_prior,
                                     effect_estimate ,
                                     list_indx_lst,
                                     low_trait,
                                     df=NULL)
{
  if(is.null(df)){
    if( is.null(G_prior$G_prior_u)){
      u_logBF <- rep(0,nrow(effect_estimate$res_f[[1]]$Bhat  ))
    }else{
      u_logBF <-  lapply(1:ncol(effect_estimate$res_u$Bhat),
                         function(k)
                           log_BFu(G_prior = G_prior$G_prior_u[[k]],
                                   Bhat    =  effect_estimate$res_u$Bhat[,k] ,
                                   Shat    =  effect_estimate$res_u$Shat[,k],
                                   low_u   =  ifelse(k%in%low_trait$low_u, TRUE,FALSE)
                           )
      )
      u_logBF <- apply(do.call(rbind, u_logBF),2,sum)
    }
    if(is.null(G_prior$G_prior_f)){
      f_logBF <- rep(0,nrow(effect_estimate$res_u[[1]] ))
    }else{
      f_logBF <- lapply( 1: length(G_prior$G_prior_f) ,function( k)
        fsusieR::log_BF(G_prior  = G_prior$G_prior_f[[k]],
                            Bhat     = effect_estimate$res_f[[k]]$Bhat,
                            Shat     = effect_estimate$res_f[[k]]$Shat,
                            indx_lst = list_indx_lst[[k]],
                            lowc_wc  = low_trait$lowc_wc[[k]]
        )
      )
      f_logBF <- apply(do.call(rbind, f_logBF),2,sum)
    }

  }else{
    ### Work here ------
    if( is.null(G_prior$G_prior_u)){
      u_logBF <- rep(0,nrow(effect_estimate$res_f[[1]]$Bhat  ))
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
      u_logBF <- apply(do.call(rbind, u_logBF),2,sum)
    }
    if(is.null(G_prior$G_prior_f)){
      f_logBF <- rep(0,nrow(effect_estimate$res_u[[1]] ))
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
      f_logBF <- apply(do.call(rbind, f_logBF),2,sum)
    }



  }

  out <- f_logBF+u_logBF
  return(out)

}


# @title Compute Log-Bayes Factor for multivariate wavelet regression
#
# @description  Compute Log-Bayes Factor for multivariate wavelet regression using a scale specific mash prior
#
# @param G_prior a scale specific mash prior
#
# @param tens_marg a list of tensor of marginal association generated by \code{cal_Bhat_Shat_tensor}
#
# @param indx_list List generated by \code{\link{gen_wavelet_indx}}
#   for the given level of resolution
#
# @return The log-Bayes factor for each covariate.
#
# @export
#
log_BF_tens <- function( G_prior, tens_marg, indx_lst){


  out <-  do.call(c,
                  lapply( 1:dim( tens_marg$tens_Bhat)[1], #number of covariate
                       function(j) cal_lbf_cov_mvsusie (G_prior,
                                                        tens_marg,
                                                        s,
                                                        j,
                                                        indx_lst
                                                        )
                       )
               )

 return(out)
}






# @title Compute posterior quantities of tensor regression for a given scale and a given covariate
#
# @description Compute posterior  quantities  of the tensor regression  using a scale specific mash prior
# @param G_prior a scale specific mash prior
#
# @param tens_marg a list of tensor of marginal association generated by \code{cal_Bhat_Shat_tensor}
# @param indx_list list generated by \code{\link{gen_wavelet_indx}} for the given level of resolution, used only with class  mixture_normal_per_scale
# @param s the scale of interest
# @param j the covariate of interest
#
# @export

get_post_effect <- function( G_prior, tens_marg, indx_lst, s,j) {
  Bhat   <-   matrix(tens_marg$tens_Bhat[j,indx_lst[[s]],],
                     ncol = dim(tens_marg$tens_Bhat)[3]
  )
  Shat   <-   matrix(tens_marg$tens_Shat[j,indx_lst[[s]],],
                     ncol = dim(tens_marg$tens_Bhat)[3]
  )
  data   <-   mash_set_data( Bhat,  Shat)

  m2     <-   mash_compute_posterior_matrices( G_prior[[s]],
                                               data
  )

  return(m2)
}





# @title Compute posterior quantities of tensor regression for a given scale
#
# @description Compute posterior  quantities  of the tensor regression  using a scale specific mash prior
# @param G_prior a scale specific mash prior
#
# @param tens_marg a list of tensor of marginal association generated by \code{cal_Bhat_Shat_tensor}
# @param indx_list list generated by \code{\link{gen_wavelet_indx}} for the given level of resolution, used only with class  mixture_normal_per_scale
# @param s the scale of interest
#
# @export
get_post_level <- function( G_prior, tens_marg, indx_lst,s, all=FALSE)
{
  res <-  lapply( 1:dim(tens_marg$tens_Bhat)[1],function(j) get_post_effect( G_prior, tens_marg, indx_lst, s,j) )

  if(!all)
  {
    PosteriorMean_level <- abind(
                                  lapply( 1:dim(tens_marg$tens_Bhat)[1],
                                          function(j) res[[j]]$PosteriorMean),
                                  along=0
                                )

    Posteriorsd_level   <- abind(
                                  lapply( 1:dim(tens_marg$tens_Bhat)[1],
                                   function(j) res[[j]]$PosteriorSD),
                                  along=0
                                )
    out <- list(PosteriorMean_level = PosteriorMean_level,
                Posteriorsd_level   = Posteriorsd_level
                )
  }else{
    PosteriorMean_level <-  abind(
                                  lapply( 1:dim(tens_marg$tens_Bhat)[1],
                                         function(j) res[[j]]$PosteriorMean
                                         ),
                                   along=0
                                  )

    Posteriorsd_level   <-  abind(
                                  lapply( 1:dim(tens_marg$tens_Bhat)[1],
                                          function(j) res[[j]]$PosteriorSD
                                          ),
                                   along=0
                                   )

    lfdr_level          <-  abind(
                                    lapply( 1:dim(tens_marg$tens_Bhat)[1],
                                           function(j) res[[j]]$lfdr
                                           ),
                                   along=0
                                   )

    NegativeProb_level  <-  abind(
                                  lapply( 1:dim(tens_marg$tens_Bhat)[1],
                                         function(j) res[[j]]$NegativeProb
                                   ),
                                   along=0)

    lfsr_level          <-  abind(lapply( 1:dim(tens_marg$tens_Bhat)[1],
                                          function(j) res[[j]]$lfsr
                                          ),
                                   along=0
                                  )

    out <- list(PosteriorMean_level = PosteriorMean_level,
                Posteriorsd_level   = Posteriorsd_level,
                lfdr_level          = lfdr_level ,
                NegativeProb_level  = NegativeProb_level,
                lfsr_level          = lfsr_level
               )
  }

  return( out)
}




# @title Compute posterior quantities of tensor regression
#
# @description Compute posterior mean and sd of the tensor regression  using a scale specific mash prior
# @param G_prior a scale specific mash prior
#
# @param tens_marg a list of tensor of marginal association generated by \code{cal_Bhat_Shat_tensor}
# @param all logical, if set as FALSE the function a=only return posterior mean and SD. If set as true return  lfdr, NegativeProb and lfsr. Set to FALSE by default
# @export
get_post_tens <- function(G_prior, tens_marg, indx_lst, all =FALSE)
{


  if ( !all){
    tt <- lapply( 1: length(indx_lst),
                  function(s)  get_post_level(G_prior, tens_marg, indx_lst,s, all =all)
    )


    post_mean_tens <- array(NA,dim =dim(tens_marg$tens_Bhat))
    post_sd_tens <- array(NA,dim =dim(tens_marg$tens_Bhat))
    for ( s in 1: 1: length(indx_lst)){
      post_mean_tens[ ,indx_lst[[s]],] <- tt[[s]]$PosteriorMean_level
      post_sd_tens  [ ,indx_lst[[s]],] <- tt[[s]]$Posteriorsd_level

    }

    out <- list(post_mean_tens = post_mean_tens,
                post_sd_tens   =   post_sd_tens
                )
  }else{


         tt <- lapply( 1: length(indx_lst),
                       function(s)  get_post_level(G_prior, tens_marg, indx_lst,s, all =all)
                      )


        post_mean_tens    <- array(NA,dim =dim(tens_marg$tens_Bhat))
        post_sd_tens      <- array(NA,dim =dim(tens_marg$tens_Bhat))
        lfdr_tens         <- array(NA,dim =dim(tens_marg$tens_Bhat))
        NegativeProb_tens <- array(NA,dim =dim(tens_marg$tens_Bhat))
        lfsr_tens         <- array(NA,dim =dim(tens_marg$tens_Bhat))


        for ( s in 1: 1: length(indx_lst)){
          post_mean_tens   [ ,indx_lst[[s]],] <- tt[[s]]$PosteriorMean_level
          post_sd_tens     [ ,indx_lst[[s]],] <- tt[[s]]$Posteriorsd_level
          lfdr_tens        [ ,indx_lst[[s]],] <- tt[[s]]$lfdr_level
          NegativeProb_tens[ ,indx_lst[[s]],] <- tt[[s]]$NegativeProb_level
          lfsr_tens        [ ,indx_lst[[s]],] <- tt[[s]]$lfsr_level

        }
        out <- list(
                   post_mean_tens    = post_mean_tens,
                   post_sd_tens      = post_sd_tens,
                   lfdr_tens         = lfdr_tens,
                   NegativeProb_tens = NegativeProb_tens,
                   lfsr_tens         = lfsr_tens
                    )
  }


 return( out)

}
#get_post_mean_tens

#get_post_sd_tens






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



# @title Compute log likelihood of multivariate gaussian model under the null
#
# @description  Compute log likelihood of multivariate gaussian model under the null
#
# @param Bhat a matrix (n_wavelet_coef x n_cond)  of MLE mean estimate
#
# @param Shat a matrix (n_wavelet_coef x n_cond)  of MLE sd estimate
#
# @return The log likelihood
#
# @export
#
mvfsusie_compute_null_loglik <- function(Bhat,Shat)
{

  out <- do.call(sum,
                 lapply( 1:nrow(Bhat),
                         function(i) dmvnorm(x=(Bhat[i,]),
                                             mean= rep(0,length(Bhat[i,])),
                                             sigma= diag((Shat[i,])^2),
                                             log=TRUE)
                 )
  )
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
             var(out$residuals)/sum(
               (X[,j]-mean(X[,j]))^2)
           )
  )
  )
}

