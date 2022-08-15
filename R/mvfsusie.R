mvfsusie <- function(Y, X, L = 2,
                  pos = NULL,

                  verbose = TRUE,
                  plot_out = TRUE,
                  maxit = 100,
                  tol = 1e-3,
                  cov_lev = 0.95,
                  min.purity=0.5,
                  lfsr_curve = 0.05,
                  filter.cs =TRUE,
                  data.driven=FALSE

)
{


  ## Input error messages

  if (is.null(pos))
  {
    pos <- 1:dim(Y[[1]])[2]
  }


  #reshaping of the data
  if ( !(length(pos)==dim(Y)[2])) #miss matching positions and number of observations
  {
    stop("Error: number of position provided different from the number of column of Y")
  }
  original_Y <-Y


  if(!is.wholenumber(log2(dim(Y)[2])) | !(sum( duplicated(diff( pos)))== (length(pos) -2)) ) #check wether dim(Y) not equal to 2^J or if the data are unevenly spaced
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

  W <- lapply(X = Y, FUN = pack_dwt)
  DW_tens <- rearrange( dwt_data, lev_res = lev_res, n_curve=dim(Y)[3])

  ### Definition of some static parameters ---
  indx_lst <-  gen_wavelet_indx(log2(length( outing_grid)))
  v1       <-  rep(1, dim(X)[1])### used in fit_lm to add a column of 1 in the design matrix

  # Wavelet transform of the inputs


  update_D <- DW_tens
  ### Definition of some dynamic parameters ---

  update_Y    <-  DW_tens #Using a column like phenotype, temporary matrix that will be regularly updated



  tens_marg  <- cal_Bhat_Shat_tensor  (DW_tens, X, v1)
  G_prior    <- init_prior_mvfsusie(tens_marg   = tens_marg,
                                    indx_list   = indx_lst,
                                    data.driven = data.driven)

  mvfsusie_obj  <- init_mvfsusie_obj (L, G_prior, DW_tens,X )

  # numerical value to check breaking condition of while
  check <- 1
  h     <- 0

  if(mvfsusie_obj$L==1)
  {
    tens_marg <- cal_Bhat_Shat_tensor  (Y, X, v1)
    tpi  <- get_pi(mvfsusie_obj,1)
    G_prior <- update_prior.mash_per_scale(G_prior, tpi= tpi ) #allow EM to start close to previous solution (to double check)
    res_EM <- EM_pi_mvfsusie(G_prior,
                             tens_marg,
                             indx_lst
    )
    susiF.obj <-  update_susiF_obj(susiF.obj = susiF.obj ,
                                   l         = 1,
                                   EM_pi     = EM_out,
                                   Bhat      = Bhat,
                                   Shat      = Shat,
                                   indx_lst  = indx_lst
    )
    susiF.obj <- update_ELBO(susiF.obj,
                             get_objective( susiF.obj = susiF.obj,
                                            Y         = Y_f,
                                            X         = X,
                                            D         = W$D,
                                            C         = W$C,
                                            indx_lst  = indx_lst
                             )
    )

  }else{
    while(check >tol & (h/L) <maxit)
    {
      for( l in 1:susiF.obj$L)
      {

        h <- h+1
        tt <- cal_Bhat_Shat(update_Y,X,v1)
        Bhat <- tt$Bhat
        Shat <- tt$Shat #UPDATE. could be nicer
        tpi <-  get_pi(susiF.obj,l)
        G_prior <- update_prior(G_prior, tpi= tpi ) #allow EM to start close to previous solution (to double check)

        EM_out  <- EM_pi(G_prior  = G_prior,
                         Bhat     =  Bhat,
                         Shat     =  Shat,
                         indx_lst =  indx_lst
        )

        susiF.obj <-  update_susiF_obj(susiF.obj = susiF.obj ,
                                       l         = l,
                                       EM_pi     = EM_out,
                                       Bhat      = Bhat,
                                       Shat      = Shat,
                                       indx_lst  = indx_lst
        )

        update_Y  <-  cal_partial_resid(
          susiF.obj = susiF.obj,
          l         = l,
          X         = X,
          D         = W$D,
          C         = W$C,
          indx_lst  = indx_lst
        )


      }#end for l in 1:L


      susiF.obj <- update_ELBO(susiF.obj,
                               get_objective( susiF.obj = susiF.obj,
                                              Y         = Y_f,
                                              X         = X,
                                              D         = W$D,
                                              C         = W$C,
                                              indx_lst  = indx_lst
                               )
      )

      sigma2    <- estimate_residual_variance(susiF.obj,Y=Y_f,X)
      susiF.obj <- update_residual_variance(susiF.obj, sigma2 = sigma2 )

      if(length(susiF.obj$ELBO)>1 )#update parameter convergence,
      {
        check <- diff(susiF.obj$ELBO)[(length( susiF.obj$ELBO )-1)]

      }
    }#end while
  }


  #preparing output
  susiF.obj <- out_prep(susiF.obj  = susiF.obj,
                        Y          = Y,
                        X          = X,
                        indx_lst   =indx_lst,
                        filter.cs  = filter.cs,
                        lfsr_curve = lfsr_curve
  )
  return(susiF.obj)
}
