

#check if one of the modality leads to some constant varaible in X
check_cst_X_sub_case <- function(X,ind_analysis){


  if( !is.null(ind_analysis$idx_f)){
  pb_f <-  do.call(c,lapply(1:length(ind_analysis$idx_f), function(k)
                          which(apply(X[ind_analysis$idx_f[[k]],],2,Rfast::colVars)==0)
      )
      )

  }else{
    pb_f<- integer(0)
  }


  if( !is.null(ind_analysis$idx_u)){
    pb_u <-  do.call(c,lapply(1:length(ind_analysis$idx_u), function(k)
      which(apply(X[ind_analysis$idx_u[[k]],],2,Rfast::colVars)==0)
    )
    )

  }else{
    pb_u<- integer(0)
  }
  return( unique (c(pb_f,pb_u)))
}


#testing if x is a wholenumber
#
#@export



is.wholenumber <- function (x, tol = .Machine$double.eps^0.5)
  abs(x - round(x)) < tol

#based on Rfast implementation
#
#@export

fast_lm <- function(x,y)
{
  be <- solve(crossprod(x),crossprod(x,y))
  resid <-  y - x %*% be
  out <- list(be = be,
              residuals = resid)
  return(out)
}


#Circular permutation on vector
# Code adapted from https://mzuer.github.io
#
#@export

shifter <- function(x, n = 1) {
  # if (n == 0) x else c(tail(x, -n), head(x, n))
  if (n == 0) x else c(utils::tail(x, n), utils::head(x, -n))
}

#shifter(c(1:10), n=-1)
# [1]  1  2  3  4  5  6  7  8  9 10
#shifter(c(1:10), n=1)
# [1] 10  1  2  3  4  5  6  7  8  9
#shifter(c(1:10), n=2)
# [1]  9 10  1  2  3  4  5  6  7  8


#@export

'%!in%' <- function(x,y)!('%in%'(x,y))



#Product bewteen a NxJ matrix and a JxKxP tensor
#returns a JxKxP tensor in which slice along dim 3 is the matrix product of the slice
#and the matrix
#@export

'%x%' <- function(mat, tens)
{
 out <-   abind::abind(
                lapply( 1:dim(tens)[3],
                        function(xi) mat%*% tens[,,xi]
                       ),
                        along =3
                )
 return(out)

}


#Product bewteen a  J vector and a JxKxP tensor
#returns a 1xKxP tensor in which slice along dim 3 is the  product between  matrix product of the slice
#and the vector

#@export

'%vxtens%' <- function(vec, tens)
{
    out <-abind::abind(
                  lapply( 1:dim(tens)[3],
                          function(xi) vec%*% tens[,,xi]
                        ),
                  along =3
                   )

  return( out)
}



#@title
#@export

# Based on Rfast implementation.
fast_lm <- function(x,y)
{

  be <- solve(crossprod(x),crossprod(x,y))
  sd <-  sqrt(Rfast::cova(y - x %*% be)/(length(x)-1))

  return(c(be,sd))
}



# @title transform 3d array into a matrix
#
# @description transform 3d array into a matrix where the number of column is equal to the length of the third dimension, code inspired from a comment of  Sven Hohenstein on stack overflow
#
# @param array  a 3 way tensor
# @return a matrix
#@export

cbind_3Darray <- function(array)
{
  #transform 3d array into a list of matrix then  concatenate each matrix finally bind them


  if(length(dim( array))==3){
    mat <- do.call(cbind, lapply ( lapply(seq(dim(array)[3]), function(x)array[ , , x]),c))
  }else{
    if(length(dim(array))==2){
      mat <- array
    }else{
      stop("Provided array is not a matrix or a 3 way tensor")
    }
  }

  return(mat)
}




# @title Check mark type for multfsusie
# @param Y list of matrices
# @param min_levres corresponds to the minimum amount of column for a trait to be considered as "functional"
# @details return of vector indicating what kind of matrices are stored in the different component of Y. USeful for multfSuSiE
is.functional <- function(Y, min_levres =4 ){


    tt2      <- c()
    dim_mark <- c()
    if( !is.null(Y$Y_f)){
      tt2 <- c(tt2,rep( 'functional', length(Y$Y_f)))
      dim_mark <- c(dim_mark, do.call( c,
                                  lapply( 1:length(Y$Y_f),
                                          function(k)
                                            ncol(Y$Y_f[[k]]))
                                  )
                    )
    }
    if( !is.null(Y$Y_u)){
      tt2 <- c(tt2, "univariate")
      dim_mark <- c(dim_mark, ncol(Y$Y_u))
    }
    ncond <- sum( ifelse( dim_mark < 2^min_levres, dim_mark, 1))
    out <- list( mark_type = tt2,
                 dim_mark  =  dim_mark,
                 ncond = ncond)


  attr(out, "class") <- 'multfsusie_data_type'
  return( out)

}

multi_array_colScale <- function(Y, scale=FALSE){

  if( !is.null(Y$Y_u))
  {
    Y$Y_u <- fsusieR:::colScale   (Y$Y_u,scale=FALSE)
  }
  if(!is.null(Y$Y_f)){
    Y$Y_f <- lapply( 1:length(Y$Y_f), function(k)  fsusieR:::colScale(Y$Y_f[[k]],scale=FALSE) )
  }

 return( Y)
}



#' @title Check mark type for multfsusie
#' @description  Allow user to define some threshold valuer for wavelet regression
#' @param thresh_u vector containing threshold for minimal variance for each univariate trait
#' @param thresh_f vector containing threshold for minimal variance for each functional trait
#' @export
threshold_set_up <- function(thresh_u, thresh_f)
{
  out <- list(thresh_u = thresh_u,
              thresh_f = thresh_f)

  return(out)
}

create_null_thresh <- function(type_mark ){

  if ( length(which(type_mark$mark_type=="univariate"))==0){#
   thresh_u <- NULL
  }else{
   thresh_u <- rep( 0, sum(type_mark$dim_mark[which(type_mark$mark_type=="univariate")]))
  }

  if ( length(which(type_mark$mark_type=="functional"))==0){#
    thresh_f <- NULL
  }else{
    thresh_f <- rep( 0, length(which(type_mark$mark_type=="functional")))
  }
 out <- threshold_set_up(thresh_u = thresh_u ,
                         thresh_f = thresh_f)
 return(out)
}

#@title function checking which
#
#@param  Y data list with two entry Y_u and Y_f containning ther differnet phenotypes
#@param thresh_lowcount an object created by \link{\code{ threshold_set_up }}

check_low_count <- function(Y, thresh_lowcount, ind_analysis ){

if(missing(ind_analysis )){
  if( !is.null(Y$Y_f)){
    temp_f <-  lapply( 1:length(Y$Y_f), function(d)
      fsusieR:::which_lowcount(Y_f=Y$Y_f[[d]] ,
                                   thresh_lowcount= thresh_lowcount$thresh_f[d]
      )
    )
  }else{
    temp_f <- NULL
  }
  if( !is.null(Y$Y_u)){
    temp_u <-  do.call(c,
                       lapply( 1:ncol(Y$Y_u), function(d)
                         (  stats::median(abs(Y$Y_u[ ,d]))<= thresh_lowcount$thresh_u[d])

                       )
    )
    temp_u <- which(temp_u)
    if(length(temp_u)==0){
      temp_u <- NULL
    }


  }else{
    temp_u <- NULL
  }
}else{


  if( !is.null(Y$Y_f)){

   temp_f <-  lapply( 1:length(Y$Y_f), function(d)
                                     fsusieR:::which_lowcount(Y_f=Y$Y_f[[d]][ind_analysis$idx_f[[d]],],
                                                           thresh_lowcount= thresh_lowcount$thresh_f[d]
                                                           )
                  )
  }else{
    temp_f <- NULL
  }
  if( !is.null(Y$Y_u)){
    temp_u <-  do.call(c,
                       lapply( 1:ncol(Y$Y_u), function(d)
                                             (  stats::median(abs(Y$Y_u[ind_analysis$idx_u[[d]],d]))<= thresh_lowcount$thresh_u[d])

                              )
                      )
    temp_u <- which(temp_u)
    if(length(temp_u)==0){
      temp_u <- NULL
    }


  }else{
    temp_u <- NULL
  }
}
  out <- list( low_wc =temp_f,
               low_u  =  temp_u)
  return( out)
}



which_notNA_pos <-  function( Y){

  if( !is.null(Y$Y_f)){
    temp_f <-  lapply( 1:length(Y$Y_f), function(d)
      which(complete.cases(Y$Y_f[[d]]))
    )

  }else{
    temp_f <- NULL
  }
  if( !is.null(Y$Y_u)){
    temp_u <-    lapply( 1:ncol(Y$Y_u), function(d)
      which(complete.cases(Y$Y_u[,d]))

    )

    if(length(temp_u)==0){
      temp_u <- NULL
    }


  }else{
    temp_u <- NULL
  }
  out <- list( idx_f =temp_f,
               idx_u  =  temp_u)
  return( out)
}


init_var_multf <- function(Y){
  sigma2          <- list()

  if(!is.null(Y$Y_f)){

    sigma2$sd_f     <- sapply(1:length(Y$Y_f) ,
                        function( k) mean(apply( Y$Y_f[[k]],
                                                 2,
                                                 function(x) stats::var(x,
                                                                        na.rm=TRUE)
                                                 )
                                          )
                             )
  }else {

    sigma2$sd_f   <- NULL
  }
  if(!is.null(Y$Y_u)){
    sigma2$sd_u     <-   apply( Y$Y_u , 2,function(x) stats::var(x,na.rm=TRUE) )


  }else{

    sigma2$sd_u   <- NULL
  }
  return(sigma2)
}
