fast_lm <- function(x,y)
{
  be <- solve(crossprod(x),crossprod(x,y))
  resid <-  y - x %*% be
  out <- list(be = be,
              residuals = resid)
  return(out)
}



#' @title transform 3d array into a matrix
#'
#' @description transform 3d array into a matrix where the number of column is equal to the length of the third dimension, code inspired from a comment of  Sven Hohenstein on stack overflow
#'
#' @param array  a 3 way tensor
#' @return a matrix
#'
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


