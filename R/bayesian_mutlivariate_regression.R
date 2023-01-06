#@title multivariate  regression
#@param Y a N by n_curve matrix
#@param x a matrix in which first column is a vector of 1 and the second n is the covariate of interest
#@param V covariance matrix of the noise, set as I if missing
#
mv_reg <- function(Y,x,V,v1)
{
  if(missing(V))#TODO correct V without estimation
  {
    V= diag(rep(1, ncol(Y)))
  }

  Bhat <-  t(solve(crossprod(x),crossprod(x,Y))) #first column contains intercept
  #second colum effect
  resid <-  Y - x %*% t(Bhat)
  out <- list( Bhat =  Bhat[,2],
               S =  var( resid)/sum(
                 (x[,2]-mean(x[,2]))^2))
  return(out)
}


#@title Bayesian multivariate  regression
#@param Y a N by n_curve matrix
#@param x a matrix in which first column is a vector of 1 and the second n is the covariate of interest
#@param V covariance matrix of the noise, set as I if missing
#@param U prior
#@importFrom mvtnorm dmvnorm
bmv_reg <- function(Y,x,V,U)
{
  if(missing(V)) #TODO correct V without estimation
  {
    V = diag(rep(1, ncol(Y)))
  }
  if(missing(U))
  {
    U =  diag(rep( sqrt(nrow(Y)) ), ncol(Y))
  }

  temp <- mv_reg(Y,x,V)

  post_Sigma1 = solve(solve(U)+ solve(temp$S))
  post_Bhat   = post_Sigma1%*%solve(temp$S)%*%temp$Bhat
  lbf         = dmvnorm(x =c(temp$Bhat) ,sigma = temp$S + U,log = TRUE) -
                     dmvnorm(x = c(temp$Bhat),sigma = temp$S,log = TRUE)

  out <- list( post_Bhat   = post_Bhat,
               post_Sigma1 = post_Sigma1,
               lbf         = lbf
              )

  return(out)
}


#@title Bayes factor for Bayesian multivariate  regression
#@param Bhat a vector of estimated effect
#@param S covariance of the estimated effect
#@param U prior
#@importFrom mvtnorm dmvnorm
multivariate_lbf = function (Bhat, S, U) {
  lbf = sapply(1:length(S),
               function(j) dmvnorm(x = Bhat[j,],sigma = S[[j]] + U,log = TRUE) -
                 dmvnorm(x = Bhat[j,],sigma = S[[j]],log = TRUE))
  lbf[which(is.nan(lbf))] = 0
  return(lbf)
}


###Checker la
