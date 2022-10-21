library(susiF.alpha)
library(mvf.susie.alpha)
library(sim1000G)

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



genotypes = SIM$gt1[1:100,] + SIM$gt2[1:100,]

(dim(genotypes))

#res <- list()
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


Y_f1 <-  matrix(rnorm((2^list_lev_res[[1]])*100 ,sd=1), nrow = 100)
Y_f2 <-  matrix(rnorm((2^list_lev_res[[2]])*100 ,sd=1), nrow = 100)

Y_u <- matrix(rnorm((n_univ)*100 ,sd=1), nrow = 100)


  tt <- sample(1:7,1)
  G <- genotypes#(tt*101):((tt+3)*102)]

  if( length(which(apply(G,2,var)==0))>0){
    G <- G[,-which(apply(G,2,var)==0)]
  }
  # G <- matrix( rnorm(100*300), nrow = 100)
  true_pos <- sample( 1:ncol(G), L)

  for ( i in 1:100){
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
  save(res, file="check_L_accuracy.RData")
}



library(susiF.alpha)
library(ashr)

tl <- list()
for (k in  1:100)
{
  #Example using curves simulated under the Mixture normal per scale prior
  rsnr <- 1 #expected root signal noise ratio
  N <- 100    #Number of individuals
  P <- 20     #Number of covariates
  pos1 <- 1   #Position of the causal covariate for effect 1
  pos2 <- 5   #Position of the causal covariate for effect 1
  lev_res <- 4#length of the molecular phenotype (2^lev_res)
  f1 <-  simu_IBSS_per_level(4 )$sim_func#first effect
  f2 <- simu_IBSS_per_level(6 )$sim_func #second effect



  G = matrix(sample(c(0, 1,2), size=N*P, replace=TRUE), nrow=N, ncol=P) #Genotype
  beta0       <- 0
  beta1       <- 1
  beta2       <- 1
  noisy.data  <- list()

  for ( i in 1:N)
  {
    f1_obs <- f1
    f2_obs <- f2
    noise <- rnorm(length(f1), sd=  2)
    noisy.data [[i]] <-   beta1*G[i,pos1]*f1_obs + beta2*G[i,pos2]*f2_obs + noise

  }
  noisy.data <- do.call(rbind, noisy.data)





  Y <- noisy.data
  X <- G
  #Running fSuSiE

  out <- susiF(Y,X,L=2 , prior = 'mixture_normal_per_scale',
               init_pi0_w = 0.5,
               control_mixsqp =  list(
                 eps = 1e-3,
                 numiter.em = 100,
                 verbose = FALSE
               )
  )
  tl[[k]]<-  out$ELBO
  print(out$sigma2)
}
for( i in 1:length(tl))
{
  plot(tl[[i]][-1])
  abline(a=max(tl[[i]]),b=0)

}
plot(out$ELBO)
diff(out$ELBO)
