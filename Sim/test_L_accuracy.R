library(susiF.alpha)
library(mvf.susie.alpha)
library(sim1000G)
set.seed(1)
examples_dir = system.file("examples", package = "sim1000G")
vcf_file = file.path(examples_dir, "region.vcf.gz")

vcf = readVCF( vcf_file , maxNumberOfVariants = 200 , min_maf = 0.02 , max_maf = NA )

# downloadGeneticMap( 4 )
readGeneticMap( chromosome = 4 )

startSimulation( vcf )



'%!in%' <- function(x,y)!('%in%'(x,y))
id = c()
for(i in 1:100) id[i] = SIM$addUnrelatedIndividual()

# Show haplotype 1  of first 5 individuals
#print(SIM$gt1[1:5,1:6])

# Show haplotype 2
#print(SIM$gt1[1:5,1:6])

N <- 30

genotypes = SIM$gt1[1:N,] + SIM$gt2[1:N,]

(dim(genotypes))

 res <- list()
#load("check_L_accuracy.RData")
for (o  in (length(res)+1):1000) {
  L <- sample(1:10, size=1)
  print(L)
  list_lev_res <- list(5,6)
  n_univ <- 3
  eff <-  list()
  for(l in 1:L){
    eff[[l]] <-   simu_effect_multfsusie (list_lev_res=list_lev_res,
                                          n_univ=n_univ, output_level = 2)
  }


Y_f1 <-  matrix(rnorm((2^list_lev_res[[1]])*N ,sd=1), nrow = N)
Y_f2 <-  matrix(rnorm((2^list_lev_res[[2]])*N ,sd=1), nrow = N)

Y_u <- matrix(rnorm((n_univ)*N ,sd=1), nrow = N)


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
  Y_f[[2]] <- Y_f1
  Y <- list( Y_f = Y_f, Y_u=Y_u)

  m1 <- multfsusie(Y=Y,
                   X=G,
                   L=11 ,
                   data.format="list_df",
                   L_start=11 ,
                   nullweight=10,
                   cal_obj =FALSE,
                   maxit=10)
  m1$cs

  true_pos[order(true_pos)]

  out <- c( length(m1$cs), #number of CS
            length(which(true_pos%in% do.call(c, m1$cs))), #number of effect found
            Reduce("+",sapply(1:length(m1$cs), function(k)
              ifelse( length(which(true_pos%in%m1$cs[[k]] ))==0, 1,0)
                          )
                   ), #number of CS without any effect
            L)
  res[[o]] <- out
  print(res)
  save(res, file="check_L_accuracy_sd1.RData")
}



 res <- list()
 #load("check_L_accuracy.RData")
 for (o  in (length(res)+1):1000) {
   L <- sample(1:10, size=1)
   print(L)
   list_lev_res <- list(5,6)
   n_univ <- 3
   eff <-  list()
   for(l in 1:L){
     eff[[l]] <-   simu_effect_multfsusie (list_lev_res=list_lev_res,
                                           n_univ=n_univ, output_level = 2)
   }


   Y_f1 <-  matrix(rnorm((2^list_lev_res[[1]])*N ,sd=2), nrow = N)
   Y_f2 <-  matrix(rnorm((2^list_lev_res[[2]])*N ,sd=2), nrow = N)

   Y_u <- matrix(rnorm((n_univ)*N ,sd=2), nrow = N)


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
   Y_f[[2]] <- Y_f1
   Y <- list( Y_f = Y_f, Y_u=Y_u)

   m1 <- multfsusie(Y=Y,
                    X=G,
                    L=11 ,
                    data.format="list_df",
                    L_start=11 ,
                    nullweight=10,
                    cal_obj =FALSE,
                    maxit=10)
   m1$cs

   true_pos[order(true_pos)]

   out <- c( length(m1$cs), #number of CS
             length(which(true_pos%in% do.call(c, m1$cs))), #number of effect found
             Reduce("+",sapply(1:length(m1$cs), function(k)
               ifelse( length(which(true_pos%in%m1$cs[[k]] ))==0, 1,0)
             )
             ), #number of CS without any effect
             L)
   res[[o]] <- out
   print(res)
   save(res, file="check_L_accuracy_sd2.RData")
 }







 res <- list()
 #load("check_L_accuracy.RData")
 for (o  in (length(res)+1):1000) {
   L <- sample(1:10, size=1)
   print(L)
   list_lev_res <- list(5,6)
   n_univ <- 3
   eff <-  list()
   for(l in 1:L){
     eff[[l]] <-   simu_effect_multfsusie (list_lev_res=list_lev_res,
                                           n_univ=n_univ, output_level = 2)
   }


   Y_f1 <-  matrix(rnorm((2^list_lev_res[[1]])*N ,sd=3), nrow = N)
   Y_f2 <-  matrix(rnorm((2^list_lev_res[[2]])*N ,sd=3), nrow = N)

   Y_u <- matrix(rnorm((n_univ)*N ,sd=3), nrow = N)


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
   Y_f[[2]] <- Y_f1
   Y <- list( Y_f = Y_f, Y_u=Y_u)

   m1 <- multfsusie(Y=Y,
                    X=G,
                    L=11 ,
                    data.format="list_df",
                    L_start=11 ,
                    nullweight=10,
                    cal_obj =FALSE,
                    maxit=10)
   m1$cs

   true_pos[order(true_pos)]

   out <- c( length(m1$cs), #number of CS
             length(which(true_pos%in% do.call(c, m1$cs))), #number of effect found
             Reduce("+",sapply(1:length(m1$cs), function(k)
               ifelse( length(which(true_pos%in%m1$cs[[k]] ))==0, 1,0)
             )
             ), #number of CS without any effect
             L)
   res[[o]] <- out
   print(res)
   save(res, file="check_L_accuracy_sd3.RData")
 }









 res <- list()
 #load("check_L_accuracy.RData")
 for (o  in (length(res)+1):1000) {
   L <- sample(1:10, size=1)
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
   Y_f[[2]] <- Y_f1
   Y <- list( Y_f = Y_f, Y_u=Y_u)

   m1 <- multfsusie(Y=Y,
                    X=G,
                    L=11 ,
                    data.format="list_df",
                    L_start=11 ,
                    nullweight=10,
                    cal_obj =FALSE,
                    maxit=10)
   m1$cs

   true_pos[order(true_pos)]

   out <- c( length(m1$cs), #number of CS
             length(which(true_pos%in% do.call(c, m1$cs))), #number of effect found
             Reduce("+",sapply(1:length(m1$cs), function(k)
               ifelse( length(which(true_pos%in%m1$cs[[k]] ))==0, 1,0)
             )
             ), #number of CS without any effect
             L)
   res[[o]] <- out
   print(res)
   save(res, file="check_L_accuracy_sd4.RData")
 }
