
#'@title Simulate function under the mash per scale prior
#'
#'@param lev_res numerical corresponds to the resolution of the simulated function (idealy between 3 and 10)
#'@param n_curve  dimension of the multivaraite time serie to generate
#'@param length_grid vector numerical corresponds to the length of the grid of sigma for mixture component(cf ash)
#'@param piO vector numerical , contain a digit  between 0 and 1, which corresponds to the null proportion ( non assocatied wavelet coefficients)
#'@param alpha numeric >0, control smoothness of the curves, should be positive and up 4 in particular d_sl ~  pi_{0,sl}  delta_0 + sum_k  pi_k N(0, 2^{- alpha * s}   sigma_k^2)
#'@param prop_decay numeric >0, control the proportion of non zero wavelet coefficient per scale, pi_{0,sl} = 1- exp(-prop_decay*s)
#'
#' @export
#' @examples
#'tt <- mvf_susie_per_level()
#'

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

    simdata = simple_sims(50,3,1)

    data = mash_set_data(simdata$Bhat, simdata$Shat)
    U.c = cov_canonical(data)

    m = mash(data, U.c, algorithm.version = 'R', posterior_samples = 100)
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
    sim_func[,j] <- wr(twav)
  }

  #plot(accessD(twav,level=6), rev( tt$D[unlist(indx_lst[7])]) )


  if(is.plot)
  {
    plot( sim_func[,1], type="l", ylim = c(min(sim_func),max(sim_func)))
    for( i in 1:n_curve)
    {
      lines(sim_func[,i], col=i)
    }
  }

  out <- list( sim_func  = sim_func,
               true_coef =tem_wt_func,
               mix_per_scale=G_level
               #emp_pi0=emp_pi0
               )
  return(out)
}


