#' @param Y list of observed time series. Length of N in which every element
#' contains a xi (number of condition) by 2^S matrix. The matrix corresponds to the
#' individuals multivariate time series
#'
#'
#' @param  all logical set to FALSE. If set to TRUE the output contains additionnal information such as
#' lfdr lfsr

mvfsusie <- function(Y, X, L = 2,
                  pos = NULL,
                  plot_out = TRUE,
                  maxit = 100,
                  tol = 1e-3,
                  cov_lev = 0.95,
                  min.purity=0.5,
                  lfsr_curve = 0.05,
                  filter.cs =TRUE,
                  data.driven=FALSE, #Still some problem with data.driven =TRUE
                  verbose= FALSE,
                  all = FALSE


)
{


  ## Input error messages

  if (is.null(pos))
  {
    pos <- 1:dim(Y[[1]])[2]
  }


  #reshaping of the data
  if ( !(length(pos)==dim(Y[[1]])[2])) #miss matching positions and number of observations
  {
    stop("Error: number of position provided different from the number of column of Y")
  }
  original_Y <-Y


  if(!is.wholenumber(log2(dim(Y[[1]])[2])) | !(sum( duplicated(diff( pos)))== (length(pos) -2)) ) #check wether dim(Y) not equal to 2^J or if the data are unevenly spaced
  {
    ##### TO REWORK ------------
    inter_pol.obj <- interpol_mat(Y, pos)
    Y             <- inter_pol.obj$Y
    bp            <- inter_pol.obj$bp
    outing_grid   <- inter_pol.obj$grid
    if(verbose)
    {
      message( "Response matrix dimensions not equal to nx 2^J \n or unevenly spaced data \n interpolation procedure used")
    }
  }  else{

    outing_grid <- 1:dim(Y[[1]])[2]
  }


  # Wavelet transform of the inputs
  W <- lapply(X = Y, FUN = pack_dwt)
  DW_tens <- rearrange( W, lev_res = lev_res, n_curve=nrow(Y[[1]]))

  ### Definition of some static parameters ---
  indx_lst <-  susiF.alpha::gen_wavelet_indx(log2(length( outing_grid)))
  v1       <-  rep(1, dim(X)[1])### used in fit_lm to add a column of 1 in the design matrix


  ### Definition of some dynamic parameters ---
  update_D <- DW_tens

  tens_marg  <- cal_Bhat_Shat_tensor  (DW_tens, X, v1)
  G_prior    <- init_prior_mvfsusie(tens_marg   = tens_marg,
                                    indx_list   = indx_lst,
                                    data.driven = data.driven,
                                    verbose     = verbose)

  mvfsusie.obj  <- init_mvfsusie_obj (L, G_prior, DW_tens,X )

  # numerical value to check breaking condition of while
  check <- 1
  h     <- 0

  if(mvfsusie.obj$L==1)
  {
    tens_marg <- cal_Bhat_Shat_tensor  (update_D , X, v1)
    tpi       <- get_pi(mvfsusie.obj,1)
    G_prior   <- update_prior.mash_per_scale(G_prior, tpi= tpi ) #allow EM to start close to previous solution (to double check)
    res_EM    <- EM_pi_mvfsusie(G_prior,
                                tens_marg,
                                indx_lst
                                )
    mvfsusie.obj <-  update_mvfsusie( mvfsusie.obj  = mvfsusie.obj ,
                                      l             = 1,
                                      EM_pi         = res_EM,
                                      tens_marg     = tens_marg,
                                      indx_lst      = indx_lst,
                                      all           = all
                                    )

  }else{
    while(check >tol & (h/L) <maxit)
    {
      for( l in 1:mvfsusie.obj$L)
      {

        tens_marg <- cal_Bhat_Shat_tensor  (Y= update_D ,X = X, v1)
        tpi       <- get_pi(mvfsusie.obj,l)
        G_prior   <- update_prior.mash_per_scale(G_prior, tpi = tpi ) #allow EM to start close to previous solution (to double check)
        res_EM    <- EM_pi_mvfsusie(G_prior,
                                    tens_marg,
                                    indx_lst
                                    )




        mvfsusie.obj <-  update_mvfsusie( mvfsusie.obj  = mvfsusie.obj ,
                                          l             = l,
                                          EM_pi         = res_EM,
                                          tens_marg     = tens_marg,
                                          indx_lst      = indx_lst,
                                          all           = all
                                        )

        #mvfsusie.obj$alpha
        update_D  <-  cal_partial_resid( mvfsusie.obj = mvfsusie.obj ,
                                         l            = l,
                                         X            = X,
                                         D            = DW_tens,
                                         indx_lst     = indx_lst
                                        )

        h <- h+1
      }#end for l in 1:L


      # mvfsusie.obj <- update_ELBO(mvfsusie.obj,
      #                         get_objective( mvfsusie.obj = mvfsusie.obj,
      #                                        Y         = Y_f,
      #                                        X         = X,
      #                                        D         = W$D,
      #                                        C         = W$C,
      #                                        indx_lst  = indx_lst
      #                         )
      #)

      sigma2       <- estimate_residual_variance(mvfsusie.obj,DW_tens,X)
      mvfsusie.obj <- update_residual_variance(mvfsusie.obj, sigma2 = sigma2 )

      if(length(mvfsusie.obj$ELBO)>1 )#update parameter convergence,
      {
        check <- diff(mvfsusie.obj$ELBO)[(length( mvfsusie.obj$ELBO )-1)]

      }

    }#end while
  }


  #preparing output
  #mvfsusie.obj <- out_prep(mvfsusie.obj = mvfsusie.obj,
  #                      Y            = Y,
  #                      X            = X,
  #                      indx_lst     = indx_lst,
  #                      filter.cs    = filter.cs
  #)
  return(mvfsusie.obj)
}
