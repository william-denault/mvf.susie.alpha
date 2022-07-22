


#noisy data list of individual multivariate time sierie
rearrange <- function( noisy.data, lev_res, n_curve)
{
  out <- c()
  for ( xi in 1:nrow(noisy.data[[1]]))
  {
    for ( ti in 1: ncol(noisy.data[[1]]))
    {
      for ( i in 1:length(noisy.data))
      {
        out <- c(out,noisy.data[[i]][xi,ti])
      }
    }

  }
  dims = c(length(noisy.data),2^lev_res, n_curve )
  return(array(out, dims))
}

#'library(ashr)
#'library(mashr)
#'library(tensorr)
#'library(wavethresh)
#'set.seed(1)
#'n_curve=3
#'lev_res=7
#'effect1 <- mvf_susie_per_level(lev_res=7,n_curve=3)$sim_func
#'effect2 <- mvf_susie_per_level(lev_res=7,n_curve=3)$sim_func
#'
#'N = 50
#'
#'#Number of covariates
#'
#'P = 10
#'
#'Choosing which variable will have an effect
#'pos1 <- 1
#'pos2 <- 2
#'
#'
#'G = matrix(sample(c(0, 1,2), size=N*P, replace=T), nrow=N, ncol=P) #Genotype
#'beta1       <- 1
#'beta2       <- 1
#'
#'
#'
#'
#'noisy.data  <- list()
#'for ( i in 1:N)
#'{
#'  f1               <- effect1
#'  f2               <- effect2
#'  noisy.data [[i]] <-  beta1*G[i,pos1]*effect1 + beta2*G[i,pos2]*effect2 + matrix( rnorm((2^lev_res)*n_curve), ncol=n_curve)
#'
#'}
#'
#'
#'noisy.data  <- list()
#'for ( i in 1:N)
#'{
#'
#'if( i >1)
#'  {beta1       <- 0
#'  beta2       <- 0
#'  }
#'  f1               <- effect1
#'  f2               <- effect2
#'  noisy.data [[i]] <-  t(beta1*G[i,pos1]*effect1 + beta2*G[i,pos2]*effect2) #+ matrix( rnorm((2^lev_res)*n_curve), ncol=n_curve)

#'}

#'tt <- array( c(noisy.data), dim = c(N, 2^lev_res, n_curve))
#'tt
#'noisy.data[[1]]
#'noisy.data[[2]]


#'temp <- rearrange( noisy.data, lev_res = 7, n_curve=3)
#'noisy.data[[1]]

#'temp[1,,1]#this corresponds to first measurement of ind 1

#'plot( temp[1,,1], noisy.data[[1]][1,])#this corresponds to first measurement of ind 1
#'plot( temp[1,,2], noisy.data[[1]][2,])#this corresponds to first measurement of ind 2
