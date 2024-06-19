
#'@title Simulate effect under the multfsusie model
#'@description Simulate effect under the multfsusie model
#'@param list_lev_res list of level of resolution (can be set to NULL if non)
#'@param n_univ number of univariate trait to be analysed
#'@param effect_univ possible effect size of the SNP on the univariate phenotype (if missing the effect size are sample at random as -1 or 1)
#'@param output_level if equal to 1 (default) provide simplified output for functional effect, setting it to 2 provide more detailled output see
#'simu_IBSS_per_level in fsusieR for additional details
#'@export
#'@examples
#'
#'
#'list_lev_res <- list(5,7)
#'n_univ <- 4
#'effect <-  simu_effect_multfsusie (list_lev_res=list_lev_res,n_univ=n_univ)
#'

simu_effect_multfsusie <- function(list_lev_res=NULL, n_univ=NULL, effect_univ,
                                   output_level=1){
      if (output_level %!in% c(1,2))
          {
            stop("output_level should be between 1 and 2")
      }
  if(missing(effect_univ)){
    effect_univ <- c( -1,1)
  }
  if( output_level==1){

    if (!is.null(list_lev_res)){
      func_effect <- list()
      for ( k in 1:length(list_lev_res)){
        func_effect[[k]] <- fsusieR::simu_IBSS_per_level  ( lev_res=list_lev_res[[k]])$sim_func
      }
    }
  }
  if(output_level==2){

    if (!is.null(list_lev_res)){
      func_effect <- list()
      for ( k in 1:length(list_lev_res)){
        func_effect[[k]] <- fsusieR::simu_IBSS_per_level  ( lev_res=list_lev_res[[k]])
      }
    }
  }

  if(!is.null(n_univ)){
    univ_effect <- sample(effect_univ, size=n_univ, replace=TRUE)
  }
  out <- list(func_effect =func_effect,
              univ_effect = univ_effect)
  return( out)
}






#'@title Simulate function under the mash per scale prior
#'@description Simulate function under the mash per scale prior
#'@param lev_res numerical corresponds to the resolution of the simulated function (idealy between 3 and 10)
#'@param n_curve  dimension of the multivaraite time serie to generate
#'@param length_grid vector numerical corresponds to the length of the grid of sigma for mixture component(cf ash)
#'@param pi0 vector numerical , contain a digit  between 0 and 1, which corresponds to the null proportion ( non assocatied wavelet coefficients)
#'@param alpha numeric >0, control smoothness of the curves, should be positive and up 4 in particular d_sl ~  pi_{0,sl}  delta_0 + sum_k  pi_k N(0, 2^{- alpha * s}   sigma_k^2)
#'@param prop_decay numeric >0, control the proportion of non zero wavelet coefficient per scale, pi_{0,sl} = 1- exp(-prop_decay*s)
#'@param is.plot logical, if true plot the simualted effect
#' @export


#to do : add decay per scale to get same level of smoothing


mvf_susie_per_level  <-function( lev_res=7,
                                 n_curve=3,
                                 length_grid= 10,
                                 pi0, #point less for the moment
                                 alpha=0.8,
                                 prop_decay=0.1,
                                 is.plot=TRUE)
{
  if(missing(pi0))
  {
    pi0 <- 1-exp(- (  prop_decay*(1:lev_res)))
  }
  #usefull to parse the indexes of wavelt coefficients

  indx_lst <- list()
  indx_lst[[1]] <- 1 #coefficient
  for ( s in 1:(lev_res-1))
  {

    indx  <- (2^s):(2^((s+1))-1)

    indx_lst[[s+1]] <- indx
  }
  #for the sake of ease same grid for each scale
  G_level <-list()
  tem_func <- matrix(0,ncol= n_curve, nrow= 2^lev_res)
  twav <- wd(tem_func[,1])
  tem_wt_func <- matrix( twav$D , length(twav$D) , n_curve )
  G_level <-list() #list of the mixture per scale
  for ( i in 1:lev_res)
  {

    simdata = mashr::simple_sims(50,3,1)

    data = mashr::mash_set_data(simdata$Bhat, simdata$Shat)
    U.c = mashr::cov_canonical(data)

    m = mashr::mash(data, U.c, algorithm.version = 'R', posterior_samples = 100)
    G_level[[i]] <- m
    for ( k in unlist(indx_lst[i]))
    {

      # clust <- sample (1:length_grid, prob=G_level[[i]]$pi, size=1)
      tem_wt_func[k,] <-   simdata$Bhat[k,]
      #print(twav$D[k])
    }


  }

  #writing the the final function
  sim_func <-  tem_func

  for( j in 1:n_curve)
  {
    twav$D <-  tem_wt_func[,j]
    sim_func[,j] <- wavethresh::wr(twav)
  }

  #plot(accessD(twav,level=6), rev( tt$D[unlist(indx_lst[7])]) )


  if(is.plot)
  {
    plot( sim_func[,1], type="l", ylim = c(min(sim_func),max(sim_func)))
    for( i in 1:n_curve)
    {
      graphics::lines(sim_func[,i], col=i)
    }
  }

  out <- list( sim_func  = sim_func,
               true_coef =tem_wt_func,
               mix_per_scale=G_level
               #emp_pi0=emp_pi0
  )
  return(out)
}


