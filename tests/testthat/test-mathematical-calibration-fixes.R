test_that("expected residual sum of squares uses X-weighted per-effect uncertainty", {
  X <- matrix(c(
    1,  0,
    0,  2,
    1, -1,
    2,  1
  ), ncol = 2, byrow = TRUE)

  alpha <- list(c(0.7, 0.3), c(0.4, 0.6))
  fitted_u <- list(
    matrix(c(0.5, -0.2), ncol = 1),
    matrix(c(-0.1, 0.3), ncol = 1)
  )
  fitted_u2 <- list(
    matrix(c(0.10, 0.05), ncol = 1),
    matrix(c(0.02, 0.08), ncol = 1)
  )
  fitted_wc <- list(
    list(matrix(c(0.3, -0.1, 0.2, 0.4), nrow = 2)),
    list(matrix(c(-0.2, 0.5, 0.1, -0.3), nrow = 2))
  )
  fitted_wc2 <- list(
    list(matrix(c(0.04, 0.03, 0.02, 0.01), nrow = 2)),
    list(matrix(c(0.02, 0.06, 0.03, 0.05), nrow = 2))
  )

  obj <- structure(list(
    L = 2L,
    alpha = alpha,
    fitted_u = fitted_u,
    fitted_u2 = fitted_u2,
    fitted_wc = fitted_wc,
    fitted_wc2 = fitted_wc2,
    n_wac = list(2L)
  ), class = "multfsusie")

  Y <- list(
    Y_u = matrix(c(0.2, -0.3, 0.5, 0.7), ncol = 1),
    Y_f = list(matrix(c(
      0.1,  0.4,
     -0.2,  0.3,
      0.6, -0.1,
      0.2,  0.8
    ), ncol = 2, byrow = TRUE))
  )

  manual_er2 <- function(y, x, mu, variance, alpha) {
    mean_by_effect <- Map(function(a, m) a * m, alpha, mu)
    second_by_effect <- Map(
      function(a, m, v) a * (m^2 + v),
      alpha, mu, variance
    )
    posterior_mean <- Reduce("+", mean_by_effect)
    d <- colSums(x^2)
    correction <- sum(vapply(seq_along(mean_by_effect), function(l) {
      second_l <- second_by_effect[[l]]
      weighted_second <- if (is.null(dim(second_l))) {
        sum(d * second_l)
      } else {
        sum(d * rowSums(second_l))
      }
      weighted_second - sum((x %*% mean_by_effect[[l]])^2)
    }, numeric(1)))
    sum((y - x %*% posterior_mean)^2) + correction
  }

  expected_u <- manual_er2(Y$Y_u[, 1], X,
                           lapply(fitted_u, function(x) x[, 1]),
                           lapply(fitted_u2, function(x) x[, 1]), alpha)
  expected_f <- manual_er2(Y$Y_f[[1]], X,
                           lapply(fitted_wc, `[[`, 1),
                           lapply(fitted_wc2, `[[`, 1), alpha)

  observed <- get_ER2(obj, Y, X)
  expect_equal(observed$uni, expected_u, tolerance = 1e-12)
  expect_equal(observed$f, expected_f, tolerance = 1e-12)

  ind_analysis <- list(idx_u = list(c(1L, 3L, 4L)),
                       idx_f = list(c(1L, 2L, 4L)))
  observed_missing <- get_ER2(obj, Y, X, ind_analysis = ind_analysis)
  expect_equal(
    observed_missing$uni,
    manual_er2(Y$Y_u[ind_analysis$idx_u[[1]], 1],
               X[ind_analysis$idx_u[[1]], , drop = FALSE],
               lapply(fitted_u, function(x) x[, 1]),
               lapply(fitted_u2, function(x) x[, 1]), alpha),
    tolerance = 1e-12
  )
  expect_equal(
    observed_missing$f,
    manual_er2(Y$Y_f[[1]][ind_analysis$idx_f[[1]], , drop = FALSE],
               X[ind_analysis$idx_f[[1]], , drop = FALSE],
               lapply(fitted_wc, `[[`, 1),
               lapply(fitted_wc2, `[[`, 1), alpha),
    tolerance = 1e-12
  )
})

test_that("quantile transformation is deterministic and preserves tied values", {
  set.seed(42)
  rng_before <- .Random.seed
  x <- c(0, 0, 0, 1, NA_real_)

  transformed <- mvf.susie.alpha:::Quantile_transform(x)

  expect_identical(.Random.seed, rng_before)
  expect_equal(transformed[1:3], rep(transformed[1], 3))
  expect_true(is.na(transformed[5]))
  expect_equal(
    transformed,
    mvf.susie.alpha:::Quantile_transform(x)
  )
})

test_that("post-hoc configurations integrate over the shared covariate", {
  logBF_trait_snp <- rbind(
    c(log(4), log(1)),
    c(log(1), log(3))
  )
  obj <- list(
    alpha = list(c(0.9, 0.1)),
    fitted_wc = NULL,
    fitted_u = list(matrix(0, nrow = 2, ncol = 2)),
    lBF_per_trait = list(list(u_logBF = logBF_trait_snp))
  )

  result <- mvf.susie.alpha:::posthoc_multfsusie(
    obj,
    prior_active = c(0.2, 0.3),
    variant_prior = c(0.6, 0.4)
  )[[1]]

  expected_bf <- c(1, 2.8, 1.8, 3.6)
  expected_prior <- c(0.56, 0.14, 0.24, 0.06)
  expected_probability <- expected_bf * expected_prior
  expected_probability <- expected_probability / sum(expected_probability)

  expect_equal(exp(result$configuration_logBF), expected_bf,
               tolerance = 1e-12)
  expect_equal(result$configuration_prior, expected_prior,
               tolerance = 1e-12)
  expect_equal(result$config_prob, expected_probability,
               tolerance = 1e-12)
  expect_equal(exp(result$logBF_trait), c(2.8, 1.8),
               tolerance = 1e-12)

  obj$alpha[[1]] <- rev(obj$alpha[[1]])
  result_with_different_alpha <- mvf.susie.alpha:::posthoc_multfsusie(
    obj,
    prior_active = c(0.2, 0.3),
    variant_prior = c(0.6, 0.4)
  )[[1]]
  expect_equal(result_with_different_alpha$config_prob, result$config_prob)
})

test_that("final filtering can return a coherent null result", {
  obj <- structure(list(
    L = 1L,
    P = 2L,
    alpha = list(c(0.6, 0.4)),
    lBF = list(c(0, 0)),
    lBF_per_trait = list(list(u_logBF = matrix(0, nrow = 1, ncol = 2))),
    cs = list(1L),
    est_pi = list(list(est_pi_u = list(c(1, 0)), est_pi_f = NULL)),
    fitted_u = list(matrix(0, nrow = 2, ncol = 1)),
    fitted_u2 = list(matrix(0, nrow = 2, ncol = 1)),
    fitted_wc = NULL,
    fitted_wc2 = NULL,
    G_prior = list(G_prior_u = list(structure(list(), class = "mixture_normal")),
                   G_prior_f = NULL),
    n_univ = 1L,
    n_wac = NULL,
    tol_null_prior = 0.001,
    pip = c(0.6, 0.4),
    KL = 1,
    lfsr_u = list(1)
  ), class = "multfsusie")

  expect_equal(unlist(get_pi0(obj, l = 1L)), 1)
  expect_equal(
    which_dummy_cs(obj, X = matrix(rnorm(8), nrow = 4), allow_empty = TRUE),
    1L
  )

  filtered <- check_cs(
    obj,
    X = matrix(rnorm(8), nrow = 4),
    allow_empty = TRUE
  )
  expect_equal(filtered$L, 0L)
  expect_length(filtered$cs, 0L)
  expect_length(filtered$alpha, 0L)
  expect_length(filtered$lBF_per_trait, 0L)
  expect_length(filtered$KL, 0L)
  expect_length(filtered$lfsr_u, 0L)
  expect_equal(filtered$pip, c(0, 0))
})
