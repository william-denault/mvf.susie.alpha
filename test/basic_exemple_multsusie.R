
library(susiF.alpha)
library(mvf.susie.alpha)
library(sim1000G)
set.seed(1)
N <-20
P <- 100
genotypes  = matrix(sample(c(0, 1,2), size=N*P, replace=T), nrow=N, ncol=P) #Genotype
res <- list()


  L <- 3
  print(L)
  list_lev_res <- list(5,6)
  n_univ <- 3
  eff <-  list()
  for(l in 1:L){
    eff[[l]] <-   simu_effect_multfsusie (list_lev_res=list_lev_res,
                                          n_univ=n_univ, output_level = 2)
  }


  Y_f1 <-  matrix(rnorm((2^list_lev_res[[1]])*N ,sd=4), nrow = N)
  Y_f2 <-  matrix(rnorm((2^list_lev_res[[2]])*N ,sd=4), nrow = N)

  Y_u <- matrix(rnorm((n_univ)*N ,sd=4), nrow = N)


  tt <- sample(1:7,1)
  G <- genotypes#(tt*101):((tt+3)*102)]

  if( length(which(apply(G,2,var)==0))>0){
    G <- G[,-which(apply(G,2,var)==0)]
  }
  # G <- matrix( rnorm(N*300), nrow = N)
  true_pos <- sample( 1:ncol(G), L)

  for ( i in 1:N){
    for ( l in 1:L){

      Y_f1[i,]<- Y_f1[i,]+eff[[l]]$func_effect[[1]]$sim_func*G[i,true_pos[[l]]]
      Y_f2[i,]<- Y_f2[i,]+eff[[l]]$func_effect[[2]]$sim_func*G[i,true_pos[[l]]]
      Y_u[i,]<- Y_u[i,]+ eff[[l]]$univ_effect*G[i,true_pos[[l]]]
    }
  }

  Y_f <- list()
  Y_f[[1]] <- Y_f1
  Y_f[[2]] <- Y_f2
  Y <- list( Y_f = Y_f, Y_u=Y_u)

  m1 <- multfsusie(Y=Y,
                   X=G,
                   L=11 ,
                   data.format="list_df",
                   L_start=11 ,
                   nullweight=10,
                   cal_obj =FALSE,
                   maxit=10, cor_small = TRUE)

  m2 <- multfsusie(Y=Y,
                   X=G,
                   L=11 ,
                   data.format="list_df",
                   L_start=11 ,
                   nullweight=10,
                   cal_obj =FALSE,
                   maxit=10 )

m1$cs

m2$cs


plot( m2$lBF[[1]], m1$lBF[[1]])


 t(m1$alpha[[1]])%*%m1$fitted_uni[[1]]
 t(m1$alpha[[2]])%*%m1$fitted_uni[[2]]
 t(m1$alpha[[3]])%*%m1$fitted_uni[[3]]
 eff[[1]]$univ_effect
 eff[[2]]$univ_effect
 eff[[3]]$univ_effect
