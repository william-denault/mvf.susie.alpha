library(mvf.susie.alpha)
library(susieRsmall)
set.seed(1)

N <- 100 # Sample size
L <- 2 # Number of effects
list_lev_res <- list(5, 6)
# Two functional phenotypes, one of length 2^5 and one of length 2^6
n_univ <- 3 # 3 univariate phenotypes
eff <- list()
for (l in 1:L) { # Simulate the multi-trait effect
  eff[[l]] <- simu_effect_multfsusie(list_lev_res=list_lev_res,
                                     n_univ=n_univ, output_level=2)
}

Y_f1 <- matrix(rnorm((length(c(eff[[l]]$func_effect[[1]]$sim_func,eff[[l]]$func_effect[[2]]$sim_func[10:30]))) * N, sd=1), nrow=N)
Y_f2 <- matrix(rnorm((length(c(eff[[l]]$func_effect[[2]]$sim_func,eff[[l]]$func_effect[[1]]$sim_func[10:30]))) * N, sd=1), nrow=N)
Y_u <- matrix(rnorm(n_univ * N, sd=1), nrow=N)
# Genotype matrix
G <- N3finemapping$X[ sample (1:500, size=N),]

G= G[ ,-which (apply(G,2,var)<0.001)]
true_pos <- sample(1:ncol(G), L) # Actually causal columns/SNPs

for (i in 1:N) {
  for (l in 1:L) {
    Y_f1[i,] <- Y_f1[i,] + c(eff[[l]]$func_effect[[1]]$sim_func,eff[[l]]$func_effect[[2]]$sim_func[10:30]) * G[i, true_pos[[l]]]
    Y_f2[i,] <- Y_f2[i,] + c(eff[[l]]$func_effect[[2]]$sim_func,eff[[l]]$func_effect[[1]]$sim_func[10:30])  * G[i, true_pos[[l]]]
    Y_u[i,]  <- Y_u[i,]  + eff[[l]]$univ_effect * G[i, true_pos[[l]]]
  }
}

Y_f <- list(Y_f1, Y_f2)
Y <- list(Y_f=Y_f, Y_u=Y_u) # Preparing data

Y1 <- list(Y_f=Y_f[[1]], Y_u=Y_u)
Y2 <- list(Y_f=Y_f[[2]], Y_u=Y_u)
pos <- list(pos1=1:ncol(Y$Y_f[[1]]), pos2=1:ncol(Y$Y_f[[2]]))
# If your signal is sampled between 1 and 64
pt= proc.time()
m1 <- multfsusie(Y=Y, X=G, pos=pos, L=3,
                 cor_small=TRUE,post_processing = "smash")
pt1= proc.time()-pt
pt= proc.time()
m11 <- fsusieR::susiF(Y=Y$Y_f[[1]], X=G, pos=pos[[1]], L=3,
                       cor_small=TRUE,
                      post_processing = "smash")
pt11= proc.time()-pt
pt= proc.time()
m12 <- fsusieR::susiF(Y=Y$Y_f[[2]], X=G, pos=pos[[2]], L=3,
                       cor_small=TRUE,
                      post_processing = "smash")
pt12= proc.time()-pt
pt= proc.time()
m21 <-   fsusieR::susiF(Y=Y$Y_f[[1]], X=G, pos=pos[[1]], L=3,
                        cor_small=FALSE,
                        post_processing = "smash")
pt21= proc.time()-pt
pt= proc.time()
m22 <-fsusieR::susiF(Y=Y$Y_f[[2]], X=G, pos=pos[[2]], L=3,
                     cor_small=FALSE,
                     post_processing = "smash")
pt22= proc.time()-pt



pt1
pt11
pt12
print(m1$cs) # Credible sets
print(m11$cs) # Credible sets

print(m12$cs) # Credible sets
print(m21$cs) # Credible sets

print(m22$cs)
print(true_pos)
