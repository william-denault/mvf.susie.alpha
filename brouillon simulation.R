library(mvf.susie.alpha)
set.seed(1234)

##############################################
# Simulation parameters
##############################################
N <- 100          # samples
P <- 50          # SNPs
n_traits <- 6     # number of functional traits
nsim_per_k <- 10  # replicates per sharing scenario
L_fit <- 3        # number of effects in fitted model
lfsr_thresh <- 0.05

# Each functional trait has a wavelet representation of resolution 6 (length 64)
list_lev_res <- rep(6, n_traits)
pos <- lapply(1:n_traits, function(i) 1:(2^list_lev_res[i]))
lev_res=6

##############################################
# Helper: Extract per-trait functional LFSR
##############################################
extract_functional_lfsr <- function(mfit) {
  # mfit$lfsr[[e]]$est_lfsr_functional is a list of length n_traits
  L <- length(mfit$lfsr)

  lfsr_mat <- matrix(NA, nrow=L, ncol=n_traits)

  for (l in 1:L) {
    lf <- mfit$lfsr[[l]]$est_lfsr_functional
    if (is.null(lf)) next
    for (t in 1:n_traits) {
      # Functional LFSR is a vector along the functional domain.
      # Collapse using min (most permissive evidence of association).
      lfsr_mat[l,t] <- min(lf[[t]], na.rm=TRUE)
    }
  }
  colnames(lfsr_mat) <- paste0("trait", 1:n_traits)
  return(lfsr_mat)
}

##############################################
# Helper: Check credible set inclusion
##############################################
is_causal_in_any_cs <- function(mfit, causal_snp) {
  if (is.null(mfit$cs)) return(NA)

  for (cs in 1:length(mfit$cs)) {

    if (causal_snp %in% mfit$cs[[cs]]) return(TRUE)
  }
  return(FALSE)
}

##############################################
# Main simulation over k = number of shared traits
##############################################
results <- list()

for (k in 1:n_traits) {
  cat("Running scenario k =", k, "shared traits\n")

  sensitivity <- numeric(nsim_per_k)
  precision   <- numeric(nsim_per_k)
  cs_cover    <- logical(nsim_per_k)

  for (rep in 1:nsim_per_k) {

    ###########################################
    # 1. Simulate genotype
    ###########################################
    G <- matrix(sample(0:2, size=N*P, replace=TRUE), nrow=N, ncol=P)
    causal_snp <- sample(1:P, 1)

    ###########################################
    # 2. Choose which functional traits are affected
    ###########################################
    affected_traits <- sort(sample(1:n_traits, k))

    ###########################################
    # 3. Simulate effects per trait
    ###########################################
    eff <- vector("list", n_traits)
    for (t in 1:n_traits) {
      if (t %in% affected_traits) {
        eff[[t]] <-    fsusieR::simu_IBSS_ash_vanilla(lev_res = lev_res)# simu_effect_fsusie()
      } else {
        eff[[t]] <- list(sim_func = rep(0, 2^lev_res))
      }
    }

    ###########################################
    # 4. Build Y_f: list of 6 functional phenotypes
    ###########################################
    Y_f <- vector("list", n_traits)
    for (t in 1:n_traits) {
      Y_f[[t]] <- matrix(rnorm(N * 2^lev_res, sd=1),
                         nrow = N, ncol = 2^lev_res)
      # Add genotype-mediated effect
      for (i in 1:N) {
        Y_f[[t]][i, ] <- Y_f[[t]][i, ] + eff[[t]]$sim_func * G[i, causal_snp]
      }
    }

    Y <- list(Y_f = Y_f, Y_u = NULL)

    ###########################################
    # 5. Fit multfsusie with HMM post-processing
    ###########################################
    mfit <-   multfsusie(Y=Y, X=G, L=L_fit, post_processing="HMM", verbose=FALSE)


    lfsr_mat <- extract_functional_lfsr(mfit)  # L x n_traits





  }

  configuration =rep(0, n_traits)
  configuration[affected_traits ]=1
  correct_SNP=is_causal_in_any_cs(mfit,causal_snp)
  results[[paste0("k",k)]] <- list(
    k=k,
    lfsr_min= lfsr_mat,
    configuration= configuration,
    correct_SNP=correct_SNP,
    causal_snp=causal_snp,
    cs_s= mfit$cs

  )
  print(results)

}





results






# Final summary table
summary_table <- do.call(rbind, lapply(results, function(r)
  data.frame(k=r$k,
             sens=r$sensitivity_mean,
             prec=r$precision_mean,
             cs=r$cs_coverage)
))

print(summary_table)

save(results, summary_table, file="multfsusie_functional_sharing_sim.RData")
