
#'@param Y a N by n_curve matrix
#'@param x a N vector
#'@param V covariance matrix of the noise, set as I if missing
mv_reg <- function(Y,x,V)
{
  if(missing(V))#TODO correct V without estimation
  {
    V= diag(rep(1, ncol(Y)))
  }
  scaling_factor <-as.numeric(crossprod(x))

  out <- list( bhat =t(Y )%*%x/scaling_factor,
               S = V/scaling_factor)
  return(out)
}



#'@param Y a N by n_curve matrix
#'@param x a N vector
#'@param V covariance matrix of the noise, set as I if missing
#'@param U prior
#'@importFrom mvtnorm dmvnorm
bmv_reg <- function(Y,x,V,U)
{
  if(missing(V)) #TODO correct V without estimation
  {
    V = diag(rep(1, ncol(Y)))
  }
  if(missing(U))
  {
    U =  diag(rep(1 ), ncol(Y))
  }

  temp <- mv_reg(Y,x,V)

  post_Sigma1 = solve(solve(U)+ solve(temp$S))
  post_bhat   = post_Sigma1%*%solve(temp$S)%*%temp$bhat
  lbf         = dmvnorm(x =c(temp$bhat) ,sigma = temp$S + U,log = TRUE) -
                     dmvnorm(x = c(temp$bhat),sigma = temp$S,log = TRUE)

  out <- list( post_bhat   = post_bhat,
               post_Sigma1 = post_Sigma1,
               lbf         = lbf
              )

  return(out)
}



multivariate_lbf = function (bhat, S, U) {
  lbf = sapply(1:length(S),
               function(j) dmvnorm(x = bhat[j,],sigma = S[[j]] + U,log = TRUE) -
                 dmvnorm(x = bhat[j,],sigma = S[[j]],log = TRUE))
  lbf[which(is.nan(lbf))] = 0
  return(lbf)
}


###Checker la
