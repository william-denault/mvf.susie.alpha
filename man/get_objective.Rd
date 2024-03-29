% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ELBO_mutlfsusie.R
\name{get_objective}
\alias{get_objective}
\alias{get_objective.multfsusie}
\title{Get objective function from data and multfsusie object}
\usage{
get_objective(multfsusie.obj, Y, X, D, C, indx_lst, ...)

\method{get_objective}{multfsusie}(multfsusie.obj, Y, X, list_indx_lst, ind_analysis, ...)
}
\arguments{
\item{multfsusie.obj}{a multfsusie object defined by \code{\link{init_multfsusie_obj}} function}

\item{Y}{Matrix of outcomes}

\item{X}{matrix of covariate}

\item{D}{matrix of wavelet D coefficients from the original input data (Y)}

\item{C}{vector of wavelet scaling coefficient from the original input data (Y)}

\item{indx_lst}{list generated by gen_wavelet_indx for the given level of resolution}
}
\value{
objective function value
}
\description{
Get objective function from data and multfsusie object
}
\keyword{internal}
