#' @param Y list of observed time series. Length of N in which every element
#' contains a xi (number of condition) by 2^S matrix. The matrix corresponds to the
#' individuals multivariate time series
#' @param X matrix of size n by p contains the covariates
#'
#' @param L the number of effect to fit (if not specified set to =2)
#'
#' @param pos vector of length J, corresponding to position/time pf
#' the observed column in Y, if missing suppose that the observation
#' are evenly spaced
#'
#' @param prior specify the prior used in susif. Three choice are
#' available "normal", "mixture_normal", "mixture_normal_per_scale"
#'
#' @param verbose If \code{verbose = TRUE}, the algorithm's progress,
#' and a summary of the optimization settings, are printed to the
#' console.
#'
#' @param plot_out If \code{plot_out = TRUE}, the algorithm's progress,
#' and a summary of the optimization settings, are ploted.
#'
#' @param tol A small, non-negative number specifying the convergence
#' tolerance for the IBSS fitting procedure. The fitting procedure
#' will halt when the difference in the variational lower bound, or
#' \dQuote{ELBO} (the objective function to be maximized), is less
#' than \code{tol}.
#'
#' @param maxit Maximum number of IBSS iterations to perform.
#'
#' @param cov_lev numeric between 0 and 1, corresponding to the
#' expected level of coverage of the cs if not specified set to 0.95
#'
#' @param min.purity minimum purity for estimated credible sets
#' @param filter.cs logical, if TRUE filter the credible set (removing low purity cs and cs with estimated prior equal to 0)
