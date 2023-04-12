library(susiF.alpha)
library(mvf.susie.alpha)
set.seed(1)
N=100
P=50

G = matrix(sample(c(0, 1,2), size=N*P, replace=T), nrow=N, ncol=P) #Genotype
beta1       <- 1
beta2       <- 1
L <- 4#actual number of effect
lf <-  list()
for(l in 1:L){
  lf[[l]] <- simu_IBSS_per_level(lev_res=5)$sim_func #functional effect for effect l
}


tt <- sample(0:4,1)


if( length(which(apply(G,2,var)==0))>0){
  G <- G[,-which(apply(G,2,var)==0)]
}
# G <- matrix( rnorm(100*300), nrow = 100)
true_pos <- sample( 1:ncol(G), L)

Y <- matrix(rnorm((2^5)*100 ,sd=4), nrow = 100)
for ( i in 1:100){
  for ( l in 1:L){
    Y[i,] <- Y[i,]+ lf[[l]]*G[i,true_pos[[l]]]
  }
}
Y_f <- list()
Y_f[[1]] <- Y

Y <- list( Y_f = Y_f, Y_u=NULL)

m1 <- multfsusie(Y=Y,
                 X=G,
                 L=11 ,
                 data.format="list_df",
                 L_start=11 ,
                 nullweight=10,
                 cal_obj =FALSE,
                 maxit=10)
plot(m1$fitted_func[[2]][[1]])
lines(lf[[3]])

sum(abs(m1$fitted_func[[2]][[1]]-lf[[3]]))

plot(m1$fitted_func[[1]][[1]])
lines(lf[[1]])
plot(m1$fitted_func[[3]][[1]])
lines(lf[[2]])
