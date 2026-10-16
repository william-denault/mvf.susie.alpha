test_that("mixture updates match ASH independently of the SNP cutoff", {
  # Identical SNP profiles have the same evidence at every cutoff. Retaining
  # more copies should not strengthen the data relative to the null penalty.
  b <- c(3, 0, -2, 0.2, 1, 4, -4, 0.1)
  Bhat <- matrix(rep(b, each = 40), nrow = 40)
  Shat <- matrix(1, 40, length(b))
  indx <- fsusieR::gen_wavelet_indx(3)
  g <- list(fitted_g = ashr::normalmix(c(0.8, 0.1, 0.1), 0, c(0, 1, 3)))
  control <- list(verbose = FALSE, eps = 1e-10, numiter.em = 40)

  for (prior in c("mixture_normal", "mixture_normal_per_scale")) {
    groups <- if (prior == "mixture_normal") list(seq_along(b)) else indx
    G_f <- structure(rep(list(g), length(groups)), class = prior)
    for (nullweight in c(1, 1.5, 10)) {
      expected_u <- ashr::ash(b[1], 1, g = g$fitted_g,
                              mixcompdist = "normal", nullweight = nullweight,
                              control = control, outputlevel = 0)$fitted_g$pi
      expected_f <- unlist(lapply(groups, function(j) {
        ashr::ash(b[j], rep(1, length(j)), g = g$fitted_g,
                  mixcompdist = "normal", nullweight = nullweight,
                  control = control, outputlevel = 0)$fitted_g$pi
      }))
      for (mode in c("univariate", "functional", "mixed")) {
        has_u <- mode != "functional"
        has_f <- mode != "univariate"
        G_prior <- structure(list(
          G_prior_u = if (has_u) list(g) else NULL,
          G_prior_f = if (has_f) list(G_f) else NULL
        ), class = "multfsusie_prior")
        estimates <- list(
          res_u = if (has_u) list(Bhat = Bhat[, 1, drop = FALSE],
                                  Shat = Shat[, 1, drop = FALSE]) else NULL,
          res_f = if (has_f) list(list(Bhat = Bhat, Shat = Shat)) else NULL
        )
        for (cutoff in c(1, 10, 40)) {
          fit <- EM_pi_multsusie(G_prior, estimates, list(indx),
                                 control_mixsqp = control,
                                 nullweight = nullweight,
                                 low_trait = list(low_u = NULL, low_wc = list(NULL)),
                                 max_SNP_EM = cutoff, tol_null_prior = 0)
          if (has_u) {
            expect_equal(as.numeric(unlist(fit$tpi_k$est_pi_u)), expected_u,
                         tolerance = 1e-6, info = paste(prior, nullweight, mode, cutoff))
          }
          if (has_f) {
            expect_equal(as.numeric(unlist(fit$tpi_k$est_pi_f)), expected_f,
                         tolerance = 1e-6, info = paste(prior, nullweight, mode, cutoff))
          }
        }
      }
    }
  }
})

test_that("the null penalty is counted once per effect and mixture", {
  pooled <- structure(c(0.5, 0.5), class = "pi_mixture_normal")
  per_scale <- structure(list(c(0.6, 0.4), c(0.7, 0.3)),
                          class = "pi_mixture_normal_per_scale")
  tpi <- list(est_pi_u = list(c(0.8, 0.2)), est_pi_f = list(pooled, per_scale))
  expected <- 9 * sum(log(c(0.8, 0.5, 0.6, 0.7)))
  expect_equal(.multfsusie_null_penalty(tpi, 10), expected)

  X <- matrix(c(-1, 0, 1), 3, 1)
  Y <- list(Y_u = matrix(c(-1, 0, 1), 3, 1), Y_f = NULL)
  obj <- structure(list(
    L = 2L, alpha = list(1, 1),
    fitted_u = list(matrix(0, 1, 1), matrix(0, 1, 1)),
    fitted_u2 = list(matrix(0, 1, 1), matrix(0, 1, 1)),
    fitted_wc = NULL, fitted_wc2 = NULL, n_wac = NULL,
    sigma2 = list(sd_u = 1, sd_f = NULL), KL = c(0, 0),
    est_pi = list(tpi, tpi)
  ), class = "multfsusie")
  ind <- list(idx_u = list(1:3), idx_f = NULL)
  expect_equal(get_objective(obj, Y, X, ind, nullweight = 10) -
                 get_objective(obj, Y, X, ind), 2 * expected)

  # No penalty must remain finite even if a mixture puts no mass on the null.
  tpi$est_pi_u[[1]] <- c(0, 1)
  expect_equal(.multfsusie_null_penalty(tpi, 1), 0)
})

test_that("univariate prior initialization uses the requested nullweight", {
  set.seed(124)
  X <- fsusieR::colScale(matrix(rnorm(30 * 6), 30, 6))
  Y <- list(Y_u = cbind(X[, 1] + rnorm(30), X[, 2] + rnorm(30)), Y_f = NULL)
  for (nullweight in c(1, 10)) {
    for (low_u in list(NULL, 2L)) {
      fit <- init_prior_multfsusie(Y, X, v1 = rep(1, 30),
                                   low_trait = list(low_u = low_u),
                                   control_mixsqp = list(verbose = FALSE),
                                   nullweight = nullweight, max_SNP_EM = 100)
      for (j in 1:2) {
        bhat <- fit$res$res_u$Bhat[, j]
        shat <- fit$res$res_u$Shat[, j]
        if (j %in% low_u) {
          bhat <- 0 * bhat
          shat <- shat + 1
        }
        expected <- ashr::ash(bhat, shat, mixcompdist = "normal",
                              nullweight = nullweight, outputlevel = 0)
        expect_equal(fit$G_prior$G_prior_u[[j]]$fitted_g, expected$fitted_g)
      }
    }
  }
})

test_that("multfsusie rejects nullweight values outside the ASH convention", {
  for (value in list(0, 0.1, -1, Inf, NA_real_, numeric(0), c(1, 10), "10")) {
    expect_error(multfsusie(nullweight = value),
                 "nullweight.*finite number at least 1")
  }
})

test_that("public fits use the same nullweight convention for both functional priors", {
  set.seed(125)
  X <- matrix(rnorm(32 * 6), 32, 6)
  Y <- list(Y_u = matrix(X[, 1] + rnorm(32), 32, 1),
            Y_f = list(matrix(rnorm(32 * 16), 32, 16) + X[, 1]))
  for (prior in c("mixture_normal", "mixture_normal_per_scale")) {
    args <- list(Y = Y, X = X, prior = prior, L = 1, L_start = 1,
                 verbose = FALSE, filter_cs = FALSE, posthoc = FALSE,
                 post_processing = "none", cal_obj = TRUE)
    default <- do.call(multfsusie, args)
    penalized <- do.call(multfsusie, c(args, list(nullweight = 10)))
    unpenalized <- do.call(multfsusie, c(args, list(nullweight = 1)))
    expect_equal(default$est_pi, penalized$est_pi)
    expect_equal(default$ELBO, penalized$ELBO)
    expect_true(all(is.finite(penalized$ELBO)))
    expect_gt(sum(abs(unlist(penalized$est_pi) - unlist(unpenalized$est_pi))), 1e-3)
  }
})

test_that("single-effect and IBSS fits monitor the penalized objective", {
  set.seed(126)
  X <- fsusieR::colScale(matrix(rnorm(32 * 6), 32, 6))
  Y <- list(Y_u = matrix(X[, 1] + rnorm(32), 32, 1), Y_f = NULL)
  centered <- list(Y_u = fsusieR::colScale(Y$Y_u, scale = FALSE), Y_f = NULL)
  ind <- list(idx_u = list(seq_len(32)), idx_f = NULL)
  for (L in 1:2) {
    fit <- multfsusie(Y, X, L = L, L_start = L, nullweight = 10,
                      standardize_Y = FALSE, residual_variance = 1,
                      greedy = FALSE, backfit = FALSE, cal_obj = TRUE,
                      verbose = FALSE, filter_cs = FALSE, posthoc = FALSE,
                      maxit = 3, tol = 1e-12)
    expect_equal(tail(fit$ELBO, 1),
                 get_objective(fit, centered, X, ind, nullweight = 10),
                 tolerance = 1e-10)
  }
})
