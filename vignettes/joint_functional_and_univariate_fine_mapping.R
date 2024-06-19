## ----include = FALSE----------------------------------------------------------
## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(mvf.susie.alpha)

## ---- include = FALSE---------------------------------------------------------

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 5,
  fig.height = 3,
  fig.align = "center",
  fig.cap = "&nbsp;",
  dpi = 175
  )
par(mar = c(1, 1, 1, 1))
 


## ----setup--------------------------------------------------------------------
library(mvf.susie.alpha)

## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = "center",
  fig.cap = "&nbsp;",
  dpi = 120
  )

 

## -----------------------------------------------------------------------------
knitr::opts_chunk$set(
	echo = TRUE,
	message = FALSE,
	warning = FALSE
)
library(fsusieR)
library(mvf.susie.alpha)
library(ashr)
library(susieR)
set.seed(1)
data(N3finemapping)
attach(N3finemapping)

rsnr <- 0.5 #expected root signal noise ratio

pos1 <- 250   #Position of the causal covariate for effect 1
pos2 <- 750   #Position of the causal covariate for effect 1
lev_res1 <- 5#length of the functional phenotype 1 phenotype (2^lev_res)

lev_res2 <- 6#length of the molecular phenotype (2^lev_res)

L <-  2

  list_lev_res <- list(lev_res1,lev_res2)
  n_univ <- 3
  eff <-  list()
  for(l in 1:L){
    eff[[l]] <-   simu_effect_multfsusie (list_lev_res=list_lev_res,
                                          n_univ=n_univ, output_level = 2)
  }
 
plot( eff[[1]]$func_effect[[1]]$sim_func, 
      type ="l",
      ylab="effect", 
      col="blue",
      main="effect in functional trait 1",
      ylim = c(-10,10)
      )
abline(a=0,b=0)
lines(eff[[2]]$func_effect[[1]]$sim_func, type="l", col="green")
legend(x=2,
       y=-3,
       lty = rep(1,2),
      legend= c("effect 1", "effect 2" ),
       col=c("blue","green" ))
plot( eff[[1]]$func_effect[[2]]$sim_func,
      type ="l",
      ylab="effect", 
      col="blue",
      main="effect in functional trait 2",
      ylim = c(-10,10)
      )
abline(a=0,b=0)
lines(eff[[2]]$func_effect[[2]]$sim_func, type="l", col="green")
legend(x=100,
       y=3,
      lty = rep(1,2),
      legend= c("effect 1", "effect 2" ),
       col=c("blue","green" ))

plot( (eff[[1]]$univ_effect+c(0,0,0.3)),#allow distinguishing the effect on the plot
     
      ylab="effect", 
      col="blue",
      main="effect on the differen univariate trait",
      ylim = c(-10,10),
      pch=19
      )
abline(a=0,b=0)
points(eff[[2]]$univ_effect,pch=19 , col="green")
legend(x=1 ,
       y=-3,
      pch = rep(19,2),
      legend= c("effect 1", "effect 2" ),
       col=c("blue","green" ))



par(mfrow=c(1,1))


