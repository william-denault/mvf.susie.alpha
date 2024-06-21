

#' @title Plot specific effect from susiF object
#'
#' @description  Plot specific effect from susiF object
#'
#' @param multfsusie.obj output of the susiF function
#'
#' @param effect  the index of the effect to be plotted
#'
#' @param cred.band logical, if TRUE, plot credible bands if
#'  susiF.obj fitted with wavelets regression. Set as TRUE by default.
#'
#' @param  lfsr.curve logical, if TRUE, plot estimated lfsr of
#' the effect at each base pair  if susiF.obj fitted with HMM regression.
#'  Set as TRUE by default.
#'
#' @param lfsr_thresh numeric, threshold for the lfsr curve, set to 0.0 by default.
#'
#' @param pip_only logical, if TRUE, plot only the pip values. Set as FALSE by default.
#'
#' @param size_line numeric, width of the plotted lines.
#'
#' @param size_point numeric, size of the point.
#'
#' @param pos_SNP vector, containing the base pair of the SNPs.
#'
#' @param point_shape vector, containing the shape of dots.
#'
#' @param title character
#' @param \dots Other arguments..
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_segment
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 geom_ribbon
#' @importFrom ggplot2 scale_fill_manual
#
#' @export
#' @examples
#' library(mvf.susie.alpha)
#' set.seed(1)
#'
#' N <- 100 # Sample size
#' P <- 100 # Number of SNPs
#' L <- 2 # Number of effects
#' list_lev_res <- list(5, 6)
#' # Two functional phenotypes, one of length 2^5 and one of length 2^6
#' n_univ <- 3 # 3 univariate phenotypes
#' eff <- list()
#' for (l in 1:L) { # Simulate the multi-trait effect
#'   eff[[l]] <- simu_effect_multfsusie(list_lev_res=list_lev_res,
#'                                      n_univ=n_univ, output_level=2)
#' }
#'
#' Y_f1 <- matrix(rnorm((2^list_lev_res[[1]]) * N, sd=1), nrow=N)
#' Y_f2 <- matrix(rnorm((2^list_lev_res[[2]]) * N, sd=1), nrow=N)
#' Y_u <- matrix(rnorm(n_univ * N, sd=1), nrow=N)
#' # Genotype matrix
#' G <- matrix(sample(c(0, 1, 2), size=N*P, replace=TRUE), nrow=N, ncol=P)
#'
#' true_pos <- sample(1:ncol(G), L) # Actually causal columns/SNPs
#'
#' for (i in 1:N) {
#'   for (l in 1:L) {
#'     Y_f1[i,] <- Y_f1[i,] + eff[[l]]$func_effect[[1]]$sim_func * G[i, true_pos[[l]]]
#'     Y_f2[i,] <- Y_f2[i,] + eff[[l]]$func_effect[[2]]$sim_func * G[i, true_pos[[l]]]
#'     Y_u[i,]  <- Y_u[i,]  + eff[[l]]$univ_effect * G[i, true_pos[[l]]]
#'   }
#' }
#'
#' Y_f <- list(Y_f1, Y_f2)
#' Y <- list(Y_f=Y_f, Y_u=Y_u) # Preparing data
#'
#' pos <- list(pos1=1:ncol(Y$Y_f[[1]]), pos2=1:ncol(Y$Y_f[[2]]))
#' # If your signal is sampled between 1 and 64
#'
#' m1 <- multfsusie(Y=Y, X=G, pos=pos, L=3)
#' plot_effect_multfsusiF(m1, effect=1, title="Effect 1")
#' plot_effect_multfsusiF(m1, effect=1, title="Effect 2")
plot_effect_multfsusiF <- function(  multfsusie.obj,
                                     effect=1,
                                     lfsr_thresh=.05,
                                     title="",
                                     cred.band = TRUE,
                                     lfsr.curve=TRUE,
                                     size_line=2,
                                     size_point=4,
                                     pos_SNP,
                                     pip_only=FALSE,
                                     point_shape, ...){

  if(missing(pos_SNP)){
    pos_SNP<-  1:length( multfsusie.obj$pip)
  }
  if( missing(point_shape)){
    point_shape <- rep( 19, length(pos_SNP))
  }
  color = c("black", "dodgerblue2", "green4", "#6A3D9A", "#FF7F00",
            "gold1", "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6",
            "#FDBF6F", "gray70", "khaki2", "maroon", "orchid1", "deeppink1",
            "blue1", "steelblue4", "darkturquoise", "green1", "yellow4",
            "yellow3", "darkorange4", "brown")
  L <-  multfsusie.obj$L

  n_wac <- multfsusie.obj$n_wac
  n_uni <- ifelse(is.null(multfsusie.obj$fitted_u[[1]]), 0,
                  length(multfsusie.obj$fitted_u[[1]])
  )






  y <-  multfsusie.obj$pip
  col_y <- rep(0, length(y))

  if (effect > L) {
    stop(paste("the specified effect should be lower or equal to ",
               L))
  }

  if( !is.null(n_wac)){

    if (cred.band& !is.null(multfsusie.obj$lfsr[[1]]$est_min_lfsr_functional)) {

      fun_effect <-  multfsusie.obj$fitted_func[[effect]]

      fun_plot <-  do.call( rbind,
                            lapply( 1: length(fun_effect), function(k)

                              data.frame(fun  = fun_effect[[k]],
                                         type = factor(rep(k, length(fun_effect[[k]]))),
                                         pos  = multfsusie.obj$outing_grid[[k]],
                                         upr  = multfsusie.obj$cred_band[[effect]][[k]][1,],
                                         lwr  = multfsusie.obj$cred_band[[effect]][[k]][2,]
                              )
                            ))


      P_func <-  ggplot( fun_plot, aes(x=pos,y=fun,ymin =  lwr ,ymax =  upr ))+
        geom_hline(yintercept = 0 , linewidth=size_line)+
        geom_line(linewidth=size_line,colour=color[effect+1])+
        facet_grid(type~., scales = "free") +
        geom_ribbon(  alpha = 0.3,fill=color[effect+1]) +
        xlab("postion") + ylab("Estimated effect")


    }
    else {



      if (lfsr.curve& !is.null( multfsusie.obj$lfsr[[1]]$est_lfsr_functional)){
        fun_effect <-  multfsusie.obj$fitted_func[[effect]]

        fun_plot <-  do.call( rbind,
                              lapply( 1: length(fun_effect), function(k)

                                data.frame(fun  = fun_effect[[k]]*ifelse(multfsusie.obj$lfsr[[effect]]$est_lfsr_functional[[k]]<lfsr_thresh,1,0),
                                           type = factor(rep(k, length(fun_effect[[k]]))),
                                           pos  = multfsusie.obj$outing_grid[[k]]

                                )
                              ))


        P_func <-  ggplot(fun_plot )+
          geom_hline(yintercept = 0.0  )+
          geom_line(  aes(x=pos,y=fun   ),linewidth=size_line,colour=color[effect+1])+

          facet_grid(type~., scales = "free") +
          xlab("postion") + ylab("Estimated effect")

      }else{

        fun_effect <-  multfsusie.obj$fitted_func[[effect]]

        fun_plot <-  do.call( rbind,
                              lapply( 1: length(fun_effect), function(k)

                                data.frame(fun= fun_effect[[k]],
                                           type=factor(rep(k, length(fun_effect[[k]]))),
                                           pos = multfsusie.obj$outing_grid[[k]]
                                )
                              ))


        P_func <-  ggplot( fun_plot, aes(x=pos,y=fun))+
          geom_hline(yintercept = 0 , linewidth=size_line)+
          geom_line(linewidth=size_line,colour=color[effect+1])+
          facet_grid(type~., scales = "free") +
          xlab("postion") + ylab("Estimated effect")


      }




    }
  }

  if (n_uni>0){

    df <- data.frame(beta= multfsusie.obj$fitted_u[[effect]],
                     lfsr = multfsusie.obj$lfsr[[effect]]$est_lfsr_univariate,
                     lfsr0 =  rep( 0,length(multfsusie.obj$fitted_u[[effect]])),
                     x=factor ( paste("trait", 1:length(multfsusie.obj$fitted_u[[effect]])))
    )
    P_uni <-  ggplot(    df ,)+
      geom_hline(yintercept = 0.05)+
      geom_point( aes(x=x,y=beta, size=size_point))+
      geom_segment(aes(x=x, xend=x, y=lfsr ,yend=lfsr0))+

      xlab("") + ylab("Estimated effect")

  }

  if( !is.null(n_wac) & n_uni>0){
    out <- gridExtra::grid.arrange(P_func,P_uni,ncol=2,top =title)
  }else{
    if(!is.null(n_wac)){
      out <- P_func
      return(out)
    }
    if(!is.null(n_wac)){
      out <- P_uni
      return(out)
    }
  }

  return(out)

}




