

# @title  wavelet transform on a matrix of functions
#
# @description function adapted from grove R package from Ma and Soriano. Perform wavelet transform of each row of a matrix
#
# @param data matrix of size NxJ where J is a power of two
#
# @param filter.number This selects the smoothness of the wavelet that you want to use in the decomposition. By default, this is 10, the Daubechies least-asymmetric orthonormal compactly supported wavelet with 10 vanishing moments. Description from Guy Nason (wavethresh package)
#
# @param family specifies the family of wavelets that you want to use. Two popular options are "DaubExPhase" and "DaubLeAsymm" but see the help for filter.select for more possibilities.. Description from Guy Nason (wavethresh package)
#
# @return A list with the following components
#
#\item{C}{ vector length n containing the wavelet C coefficient}
#\item{D}{ matrix of size nx 2^(J -1) where each row contains the wavelet D coefficients, ordered in the same way as in the wavethresh package}
#\item{family}{ used for the wavelet transform}
#\item{filter.number}{ }
#
# @importFrom wavethresh accessC
# @importFrom wavethresh wd
# @importFrom stats complete.cases
# @export
DWT2 <- function (data, filter.number = 10, family = "DaubLeAsymm")
{

  NA_pos <- which(!stats::complete.cases(data))
  if (length(NA_pos )>0){
    data [is.na(data)]<-0
  }
  J <- ncol(data)
  n <- nrow(data)
  D <- matrix(NA, nrow = n, ncol = J - 1)
  C <- rep(NA, n)
  for (i in 1:n) { ## Speed Gain
    temp <- wavethresh::wd(data[i, ], filter.number = filter.number,
                           family = family)
    D[i, ] <- temp$D
    C[i] <- wavethresh::accessC(temp, level = 0)
  }
  if (length(NA_pos )>0){
    D [NA_pos,] <- NA
    C [NA_pos]  <- NA
  }

  output <- list(C = C, D = D, J = log2(J), filter.number = filter.number,
                 family = family)
  class(output) <- "DWT"
  return(output)
}




# @title Compact code for wavelet transform on a matrix of functions
#
# @description function adapted from grove R package from Ma and Soriano. Perform wavelet transform of each row of a matrix
#
# @param Y matrix of size NxJ where J is a power of two
#
# @return A Matrix in which the C compoent is stored in the last column
# @export

pack_dwt <- function( Y)
{
  W <- DWT2(Y)
  return(cbind( W$D,W$C))
}


