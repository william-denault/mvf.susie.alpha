

#check if one of the modality leads to some constant varaible in X
check_cst_X_sub_case <- function(X, ind_analysis) {
  X <- as.matrix(X)
  if (ncol(X) == 0L) {
    return(integer(0))
  }

  row_groups <- c(ind_analysis$idx_f, ind_analysis$idx_u)
  if (length(row_groups) == 0L) {
    row_groups <- list(seq_len(nrow(X)))
  }

  is_constant <- rep(FALSE, ncol(X))
  for (rows in row_groups) {
    rows <- unique(as.integer(rows))
    rows <- rows[!is.na(rows) & rows >= 1L & rows <= nrow(X)]
    if (length(rows) < 2L) {
      stop("At least two analyzed observations are required for every modality")
    }
    variance <- apply(
      X[rows, , drop = FALSE],
      2,
      stats::var
    )
    is_constant <- is_constant | !is.finite(variance) | variance == 0
  }

  which(is_constant)
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

.normalize_standardize_Y <- function(standardize_Y, Y) {
  n_f <- if (is.null(Y$Y_f)) 0L else length(Y$Y_f)
  n_u <- if (is.null(Y$Y_u)) 0L else ncol(Y$Y_u)

  normalize_component <- function(x, n, label) {
    if (n == 0L) {
      if (!is.null(x) && length(x) > 0L) {
        stop(paste0("`standardize_Y$", label,
                    "` was supplied but that outcome block is absent"))
      }
      return(logical())
    }
    if (is.null(x)) {
      return(rep(TRUE, n))
    }
    if (!is.logical(x) || anyNA(x) || !length(x) %in% c(1L, n)) {
      stop(paste0("`standardize_Y$", label,
                  "` must be TRUE/FALSE or have one logical value per modality"))
    }
    rep_len(x, n)
  }

  if (is.logical(standardize_Y) && length(standardize_Y) == 1L &&
      !is.na(standardize_Y)) {
    return(list(Y_f = rep(standardize_Y, n_f),
                Y_u = rep(standardize_Y, n_u)))
  }
  if (!is.list(standardize_Y) || is.null(names(standardize_Y)) ||
      any(names(standardize_Y) == "") ||
      any(!names(standardize_Y) %in% c("Y_f", "Y_u")) ||
      anyDuplicated(names(standardize_Y))) {
    stop(paste0(
      "`standardize_Y` must be TRUE/FALSE or a named list with ",
      "components `Y_f` and/or `Y_u`"
    ))
  }

  list(
    Y_f = normalize_component(standardize_Y$Y_f, n_f, "Y_f"),
    Y_u = normalize_component(standardize_Y$Y_u, n_u, "Y_u")
  )
}


.normalize_residual_variance <- function(residual_variance, Y) {
  n_f <- if (is.null(Y$Y_f)) 0L else length(Y$Y_f)
  n_u <- if (is.null(Y$Y_u)) 0L else ncol(Y$Y_u)

  normalize_component <- function(x, n, label) {
    if (n == 0L) {
      if (!is.null(x) && length(x) > 0L) {
        stop(paste0("`residual_variance$", label,
                    "` was supplied but that outcome block is absent"))
      }
      return(numeric())
    }
    if (is.null(x)) {
      return(rep(NA_real_, n))
    }
    if (!is.numeric(x) || !length(x) %in% c(1L, n)) {
      stop(paste0("`residual_variance$", label,
                  "` must be numeric and have one value per modality"))
    }
    x <- rep_len(as.numeric(x), n)
    if (any(!is.na(x) & (!is.finite(x) | x <= 0))) {
      stop(paste0("Non-missing values in `residual_variance$", label,
                  "` must be finite and positive"))
    }
    x
  }

  if (is.null(residual_variance)) {
    return(list(sd_f = rep(NA_real_, n_f),
                sd_u = rep(NA_real_, n_u)))
  }
  if (is.numeric(residual_variance) && length(residual_variance) == 1L) {
    residual_variance <- list(Y_f = residual_variance,
                              Y_u = residual_variance)
    if (n_f == 0L) residual_variance$Y_f <- NULL
    if (n_u == 0L) residual_variance$Y_u <- NULL
  }
  if (!is.list(residual_variance) || is.null(names(residual_variance)) ||
      any(names(residual_variance) == "") || anyDuplicated(names(residual_variance))) {
    stop(paste0(
      "`residual_variance` must be NULL, a positive scalar, or a named list ",
      "with components `Y_f` and/or `Y_u`"
    ))
  }

  aliases <- c(Y_f = "Y_f", Y_u = "Y_u", sd_f = "Y_f", sd_u = "Y_u")
  if (any(!names(residual_variance) %in% names(aliases))) {
    stop(paste0(
      "`residual_variance` must use components `Y_f` and/or `Y_u` ",
      "(`sd_f` and `sd_u` are accepted aliases)"
    ))
  }
  canonical_names <- unname(aliases[names(residual_variance)])
  if (anyDuplicated(canonical_names)) {
    stop("Specify each residual-variance outcome block only once")
  }
  names(residual_variance) <- canonical_names

  list(
    sd_f = normalize_component(residual_variance$Y_f, n_f, "Y_f"),
    sd_u = normalize_component(residual_variance$Y_u, n_u, "Y_u")
  )
}


.validate_multfsusie_response_controls <- function(standardize_Y,
                                                    residual_variance,
                                                    quantile_trans = FALSE) {
  fixed_f <- !is.na(residual_variance$sd_f)
  fixed_u <- !is.na(residual_variance$sd_u)
  if (any(fixed_f & standardize_Y$Y_f) ||
      any(fixed_u & standardize_Y$Y_u)) {
    stop("Fixed `residual_variance` requires `standardize_Y = FALSE` for the same modality")
  }
  if (isTRUE(quantile_trans) && (any(fixed_f) || any(fixed_u))) {
    stop("Fixed `residual_variance` is incompatible with `quantile_trans`")
  }
  invisible(TRUE)
}


.apply_fixed_residual_variance <- function(sigma2, residual_variance) {
  for (component in c("sd_f", "sd_u")) {
    requested <- residual_variance[[component]]
    if (length(requested) == 0L) {
      next
    }
    if (is.null(sigma2[[component]]) || length(sigma2[[component]]) != length(requested)) {
      stop(paste0("Residual-variance component `", component,
                  "` is not aligned with the outcomes"))
    }
    fixed <- !is.na(requested)
    sigma2[[component]][fixed] <- requested[fixed]
  }
  sigma2
}


.outcome_scale <- function(x) {
  scale <- attr(x, "scaled:scale")
  if (is.null(scale)) rep(1, ncol(x)) else as.numeric(scale)
}


multi_array_colScale <- function(Y, scale=TRUE){
  standardize_Y <- .normalize_standardize_Y(scale, Y)

  if (!is.null(Y$Y_u)) {
    selected <- standardize_Y$Y_u
    if (all(selected) || !any(selected)) {
      # Keep the pre-existing all-TRUE and all-FALSE paths exactly unchanged.
      Y$Y_u <- fsusieR::colScale(Y$Y_u, scale = selected[1])
    } else {
      center <- colMeans(Y$Y_u, na.rm = TRUE)
      scale_u <- rep(1, ncol(Y$Y_u))
      scale_u[selected] <- matrixStats::colSds(
        Y$Y_u[, selected, drop = FALSE],
        center = center[selected],
        na.rm = TRUE
      )
      scale_u[!is.finite(scale_u) | scale_u == 0] <- 1
      Y$Y_u <- sweep(Y$Y_u, 2, center, "-")
      Y$Y_u <- sweep(Y$Y_u, 2, scale_u, "/")
      attr(Y$Y_u, "scaled:center") <- center
      attr(Y$Y_u, "scaled:scale") <- scale_u
    }
  }
  if (!is.null(Y$Y_f)) {
    Y$Y_f <- lapply(seq_along(Y$Y_f), function(k) {
      fsusieR::colScale(Y$Y_f[[k]], scale = standardize_Y$Y_f[k])
    })
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

Quantile_transform <- function(x, ties.method = "average")
{
  if (!ties.method %in% c("average", "random", "max", "min", "first", "last")) {
    stop("Unsupported ties.method for the quantile transformation")
  }

  observed <- !is.na(x)
  out <- rep(NA_real_, length(x))
  n_observed <- sum(observed)

  if (n_observed == 0L) {
    return(out)
  }

  x_rank <- rank(x[observed], ties.method = ties.method)
  out[observed] <- stats::qnorm((x_rank - 0.5) / n_observed)
  return(out)
}

mfsusie_Quantile_transform=function(Y){

  if (!is.null(Y$Y_f)){
    for ( k in 1: length(Y$Y_f))
      Y$Y_f[[k]]=apply(Y$Y_f[[k]],2,  Quantile_transform )
  }
  if( !is.null(Y$Y_u)){
    Y$Y_u=apply(Y$Y_u,2,  Quantile_transform )
  }

  return(Y)

}


#@title function checking which
#
#@param  Y data list with two entry Y_u and Y_f containning ther differnet phenotypes
#@param thresh_lowcount an object created by \link{\code{ threshold_set_up }}

check_low_count <- function(Y, thresh_lowcount, ind_analysis,type_mark ){



  if (length(thresh_lowcount )==1 & is.numeric(thresh_lowcount)){
    threshs <- create_null_thresh(type_mark = type_mark)
    if ( !is.null(threshs$thresh_u)){
      threshs$thresh_u= rep(thresh_lowcount, length(threshs$thresh_u))
    }
    if ( !is.null(threshs$thresh_f)){
      threshs$thresh_f= rep(thresh_lowcount, length(threshs$thresh_f))
    }
    thresh_lowcount=threshs
  }

if(missing(ind_analysis )){
  if( !is.null(Y$Y_f)){
    temp_f <-  lapply( 1:length(Y$Y_f), function(d)
      fsusieR::which_lowcount(Y_f=Y$Y_f[[d]] ,
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
                                     fsusieR::which_lowcount(Y_f=Y$Y_f[[d]][ind_analysis$idx_f[[d]],],
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


.normalize_function_positions <- function(Y_f, pos = NULL) {
  n_modality <- length(Y_f)
  if (is.null(pos)) {
    pos <- vector("list", n_modality)
  }
  if (!is.list(pos)) {
    stop("pos must be a list with one entry per functional modality")
  }
  if (length(pos) > n_modality) {
    stop("pos has more entries than Y$Y_f")
  }
  if (length(pos) < n_modality) {
    length(pos) <- n_modality
  }

  for (i in seq_along(Y_f)) {
    if (is.null(pos[[i]])) {
      pos[[i]] <- seq_len(ncol(Y_f[[i]]))
    }
    if (length(pos[[i]]) != ncol(Y_f[[i]])) {
      stop(paste(
        "Error: number of position provided different from the number of column of Y$Y_f, entry",
        i
      ))
    }
    if (!is.numeric(pos[[i]]) || anyNA(pos[[i]]) ||
        any(!is.finite(pos[[i]])) || any(diff(pos[[i]]) <= 0)) {
      stop(paste(
        "Error: positions must be finite and strictly increasing for Y$Y_f, entry",
        i
      ))
    }
  }

  pos
}


init_var_multf <- function(Y){
  sigma2          <- list()
  variance_floor <- sqrt(.Machine$double.eps)

  if(!is.null(Y$Y_f)){

    sigma2$sd_f     <- pmax(sapply(1:length(Y$Y_f) ,
                        function( k) mean(apply( Y$Y_f[[k]],
                                                 2,
                                                 function(x) stats::var(x,
                                                                        na.rm=TRUE)
                                                 )
                                          )
                             ), variance_floor)
  }else {

    sigma2$sd_f   <- NULL
  }
  if(!is.null(Y$Y_u)){
    sigma2$sd_u     <- pmax(
      apply(Y$Y_u, 2, function(x) stats::var(x, na.rm=TRUE)),
      variance_floor
    )


  }else{

    sigma2$sd_u   <- NULL
  }
  return(sigma2)
}




.log_sum_exp <- function(x) {
  max_x <- max(x)
  if (!is.finite(max_x)) {
    return(max_x)
  }
  max_x + log(sum(exp(x - max_x)))
}

.normalize_variant_prior <- function(variant_prior, n_variant) {
  if (is.null(variant_prior)) {
    return(rep(1 / n_variant, n_variant))
  }
  if (!is.numeric(variant_prior) || length(variant_prior) != n_variant ||
      anyNA(variant_prior) || any(!is.finite(variant_prior)) ||
      any(variant_prior < 0) ||
      sum(variant_prior) <= 0) {
    stop("variant_prior must be a non-negative numeric vector with one value per covariate")
  }
  variant_prior / sum(variant_prior)
}

get_cs_logBF_multfsusie <- function(alpha_l = NULL, logBF_trait_snp,
                                    variant_prior = NULL) {
  logBF_trait_snp <- as.matrix(logBF_trait_snp)
  variant_prior <- .normalize_variant_prior(variant_prior,
                                            ncol(logBF_trait_snp))
  log_variant_prior <- log(variant_prior)

  apply(logBF_trait_snp, 1L,
        function(x) .log_sum_exp(log_variant_prior + x))
}
