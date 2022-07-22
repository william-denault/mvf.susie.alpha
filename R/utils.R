fast_lm <- function(x,y)
{
  be <- solve(crossprod(x),crossprod(x,y))
  resid <-  y - x %*% be
  out <- list(be = be,
              residuals = resid)
  return(out)
}



#' @title concatenate 3d array into a
#'
#' @description concatenate 3d array into a matrix, code inspired from a comment of  Sven Hohenstein on stack overflow
#'
#' @param array  a 3 way tensor
#' @return a matrix
#'
rbind_3Darray <- function(array)
{
  #transform 3d array into a list of matrix then  concatenate the list
  if(length(dim( array))==3){
    mat <- do.call(rbind, lapply(seq(dim(array)[3]), function(x)array[ , , x]))
  }else{
    if(length(dim(array))==2){
      mat <- array
    }else{
      stop("Provided array is not a matrix or a 3 way tensor")
    }
  }

  return(mat)
}
