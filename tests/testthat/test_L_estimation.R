library(testthat)
library(susieR)
library(fsusieR)

set.seed(123)

test_that("susiF correctly estimates number of effects when nullweight is high enough", {

  N <- 100  # sample size
  P <- 100  # number of SNPs
  L <- 2    # number of effects

  # Genotype matrix (assuming N3finemapping$X exists in test env)
  G <- N3finemapping$X[sample(1:nrow(N3finemapping$X), size = N), ]
  # drop quasi-constant columns
  bad <- which(apply(G, 2, var) < 1e-4)
  if (length(bad) > 0) {
    G <- G[, -bad, drop = FALSE]
  }

  # simulate multi-resolution functional effects
  list_lev_res <- list(5, 6)
  n_univ <- 3
  eff <- vector("list", L)
  for (l in 1:L) {
    eff[[l]] <- simu_effect_multfsusie(
      list_lev_res = list_lev_res,
      n_univ = n_univ,
      output_level = 2
    )
  }

  ## lengths for the two functional phenotypes
  len_f1 <- length(c(
    eff[[1]]$func_effect[[1]]$sim_func,
    eff[[1]]$func_effect[[2]]$sim_func[10:30]
  ))
  len_f2 <- length(c(
    eff[[1]]$func_effect[[2]]$sim_func,
    eff[[1]]$func_effect[[1]]$sim_func[10:30]
  ))

  # simulate noise-only responses
  Y_f1 <- matrix(rnorm(len_f1 * N, sd = 1), nrow = N)
  Y_f2 <- matrix(rnorm(len_f2 * N, sd = 1), nrow = N)
  Y_u  <- matrix(rnorm(n_univ * N, sd = 1), nrow = N)

  # choose causal SNPs
  true_pos <- sample(1:ncol(G), L)

  # add genetic effects to the phenotypes
  for (i in 1:N) {
    for (l in 1:L) {
      Y_f1[i, ] <- Y_f1[i, ] +
        c(eff[[l]]$func_effect[[1]]$sim_func,
          eff[[l]]$func_effect[[2]]$sim_func[10:30]) *
        G[i, true_pos[l]]

      Y_f2[i, ] <- Y_f2[i, ] +
        c(eff[[l]]$func_effect[[2]]$sim_func,
          eff[[l]]$func_effect[[1]]$sim_func[10:30]) *
        G[i, true_pos[l]]

      Y_u[i, ]  <- Y_u[i, ] + eff[[l]]$univ_effect * G[i, true_pos[l]]
    }
  }

  # run fsusie separately on the two functional traits
  m1 <- susiF(
    Y = Y_f1,
    X = G,
    L = 3,
    verbose = FALSE,
    control_mixsqp = list(
      verbose = FALSE,
      eps = 1e-6,
      numiter.em = 40
    ),
    nullweight = .1
  )

  m2 <- susiF(
    Y = Y_f2,
    X = G,
    L = 3,
    verbose = FALSE,
    control_mixsqp = list(
      verbose = FALSE,
      eps = 1e-6,
      numiter.em = 40
    ),
    nullweight = .1
  )

  # count credible sets that do NOT contain any true causal SNP
  n_false_fsusie1 <- Reduce(
    "+",
    sapply(seq_along(m1$cs), function(k) {
      ifelse(length(which(true_pos %in% m1$cs[[k]])) == 0, 1, 0)
    })
  )

  n_false_fsusie2 <- Reduce(
    "+",
    sapply(seq_along(m2$cs), function(k) {
      ifelse(length(which(true_pos %in% m2$cs[[k]])) == 0, 1, 0)
    })
  )

  # tests
  expect_equal(n_false_fsusie1, 0)
  expect_lte(n_false_fsusie2, 0)  # i.e. also 0
})
