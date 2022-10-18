#'@param list_lev_res list of level of resoultion (can be set to NULL if non)
#'@param n_univ number of univariate trait to be analysed

simu_effec_multfsusie <- function(list_lev_res=NULL, n_univ=NULL){


  if (!is.null(list_lev_res)){
    func_effect <- list()
    for ( k in 1:length(list_lev_res)){
      func_effect[[k]] <- simu_IBSS_per_level  ( lev_res=list_lev_res[[k]])
    }
  }

  if(!is.null(n_univ)){
    univ_effect <- sample(c(-1,1), size=n_univ, replace=TRUE)
  }
  out <- list(func_effect =func_effect,
              univ_effect = univ_effect)
  return( out)
}
