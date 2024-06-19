rm(list=ls()
)
library(mvf.susie.alpha)
library(testthat)
set.seed(1)
N <- 100  #Sample size
P= 10 # number of SNP
L <-2
list_lev_res <- list(5,6) # two functional phenotypes ,
#one of length 2^5, and one of length 2^6)
n_univ <- 3 #3 univariate phenotypes
eff <-  list()
for(l in 1:L){ #Simulate the mult-trait effect
  eff[[l]] <-   simu_effect_multfsusie (list_lev_res=list_lev_res,
                                        n_univ=n_univ, output_level = 2)
}



eff[[1]]$func_effect[[1]]$sim_func <- eff[[l]]$func_effect[[1]]$sim_func*0
eff[[1]]$func_effect[[1]]$sim_func [ 10] <- 2


eff[[2]]$func_effect[[1]]$sim_func <- eff[[l]]$func_effect[[1]]$sim_func*0
eff[[2]]$func_effect[[1]]$sim_func [ 20] <- 2

Y_f1 <-  matrix(rnorm((2^list_lev_res[[1]])*N ,sd=5), nrow = N)
Y_f2 <-  matrix(rnorm((2^list_lev_res[[2]])*N ,sd=5), nrow = N)

Y_u <- matrix(rnorm((n_univ)*N ,sd=1), nrow = N)


G = matrix(sample(c(0, 1,2), size=N*P, replace=TRUE), nrow=N, ncol=P) #Genotype





true_pos <- sample( 1:ncol(G), L)# actually causal column/SNP

for ( i in 1:N){
  for ( l in 1:L){

    Y_f1[i,]<- Y_f1[i,]+eff[[l]]$func_effect[[1]]$sim_func*G[i,true_pos[[l]]]
    Y_f2[i,]<- Y_f2[i,]+eff[[l]]$func_effect[[2]]$sim_func*G[i,true_pos[[l]]]
    Y_u[i,]<- Y_u[i,]+ eff[[l]]$univ_effect*G[i,true_pos[[l]]]
  }
}

Y_f <- list()

Y_f1[1:20,1] <-NA ### heree pb because NA in the first modality
X <- G
X[,1] <- 1
X[1:20,1] <- 2### and X constant in the first modality


Y_f[[1]] <- Y_f1
Y_f[[2]] <- Y_f2
Y <- list( Y_f = Y_f, Y_u=Y_u) # preparing data ,
#current ouput type expect list of two which element named
#Y_f for functional trait and Y_u for univariate trait


ind_analysis <- which_notNA_pos(Y)
#remove column of X constant in some sub cases
tidx <- check_cst_X_sub_case(X,ind_analysis)

print(tidx)


pos = list(pos1= 1: ncol(Y$Y_f[[1]]),
           pos2= 1: ncol(Y$Y_f[[2]])) # if you signal is sample between 1 and 64

m1 <- multfsusie(Y=Y,
                 X=X,
                 pos=pos,
                 L=11 ,
                 nullweight=10,
                 maxit=10 ,
                 post_processing = "HMM")
test_that("check if dummy colum actually removed", { expect_equal(length(m1$pip), ncol(X)-1)})


