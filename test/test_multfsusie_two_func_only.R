library(susiF.alpha)
library(mvf.susie.alpha)
N=100
P=50
set.seed(2)
G = matrix(sample(c(0, 1,2), size=N*P, replace=T), nrow=N, ncol=P) #Genotype
beta1       <- 1
beta2       <- 1
L <- 4#actual number of effect
lf <-  list()
lf2 <-  list()

for(l in 1:L){
  lf [[l]] <- simu_IBSS_per_level(lev_res=5)$sim_func #functional effect for effect l
  lf2[[l]] <- simu_IBSS_per_level(lev_res=7)$sim_func #functional effect for effect l

}


tt <- sample(0:4,1)


if( length(which(apply(G,2,var)==0))>0){
  G <- G[,-which(apply(G,2,var)==0)]
}
# G <- matrix( rnorm(100*300), nrow = 100)
true_pos <- sample( 1:ncol(G), L)

Y1 <- matrix(rnorm((2^5)*100 ,sd=1), nrow = 100)
Y2 <- matrix(rnorm((2^7)*100 ,sd=1), nrow = 100)

for ( i in 1:100){
  for ( l in 1:L){
    Y1[i,] <- Y1[i,]+ lf[[l]]*G[i,true_pos[[l]]]
    Y2[i,] <- Y2[i,]+ lf2[[l]]*G[i,true_pos[[l]]]
  }
}
Y_f <- list()
Y_f[[1]] <- Y1

Y <- list( Y_f = Y_f, Y_u=NULL)

m1 <- multfsusie(Y=Y,
                 X=G,
                 L=11 ,

                 L_start=11 ,
                 nullweight=10,
                 cal_obj =FALSE,
                 maxit=10)
plot(m1$fitted_func[[1]][[1]])
lines(lf[[2]])
plot(m1$fitted_func[[2]][[1]])
lines(lf[[4]])
plot(m1$fitted_func[[3]][[1]])
lines(lf[[3]])
plot(m1$fitted_func[[4]][[1]])
lines(lf[[1]])


Y_f <- list()
Y_f[[1]] <- Y1
Y_f[[2]] <- Y2


Y <- list( Y_f = Y_f, Y_u=NULL)

m2 <- multfsusie(Y=Y,
                 X=G,
                 L=6 ,

                 L_start=6 ,
                 nullweight=10,
                 cal_obj =TRUE,
                 maxit=10)
plot(m2$fitted_func[[1]][[1]], type="l", col="red")
lines(lf[[2]])
lines(m1$fitted_func[[1]][[1]], col="green")
plot(m2$fitted_func[[2]][[1]], type="l", col="red")
lines(lf[[4]])
lines(m1$fitted_func[[2]][[1]], col="green")

plot(m2$fitted_func[[3]][[1]], type="l", col="red")
lines(lf[[3]])
lines(m1$fitted_func[[3]][[1]], col="green")
plot(m2$fitted_func[[4]][[1]], type="l", col="red")
lines(lf[[1]])
lines(m1$fitted_func[[4]][[1]], col="green")

