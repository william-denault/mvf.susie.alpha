test_that("response controls are normalized per modality", {
  Y <- list(
    Y_f = list(matrix(0, 4, 8), matrix(0, 4, 16)),
    Y_u = matrix(0, 4, 3)
  )

  expect_identical(
    mvf.susie.alpha:::.normalize_standardize_Y(TRUE, Y),
    list(Y_f = c(TRUE, TRUE), Y_u = c(TRUE, TRUE, TRUE))
  )
  standardize <- mvf.susie.alpha:::.normalize_standardize_Y(
    list(Y_f = c(FALSE, TRUE), Y_u = c(TRUE, FALSE, TRUE)),
    Y
  )
  expect_identical(
    standardize,
    list(Y_f = c(FALSE, TRUE), Y_u = c(TRUE, FALSE, TRUE))
  )

  residual_variance <- mvf.susie.alpha:::.normalize_residual_variance(
    list(Y_f = c(0.4, NA), Y_u = c(NA, 2, NA)),
    Y
  )
  expect_equal(residual_variance$sd_f, c(0.4, NA))
  expect_equal(residual_variance$sd_u, c(NA, 2, NA))
  expect_silent(
    mvf.susie.alpha:::.validate_multfsusie_response_controls(
      standardize,
      residual_variance
    )
  )

  expect_error(
    mvf.susie.alpha:::.normalize_standardize_Y(
      list(Y_f = c(TRUE, FALSE, TRUE)),
      Y
    ),
    "one logical value per modality"
  )
  expect_error(
    mvf.susie.alpha:::.validate_multfsusie_response_controls(
      list(Y_f = c(TRUE, TRUE), Y_u = c(TRUE, TRUE, TRUE)),
      residual_variance
    ),
    "same modality"
  )
  expect_error(
    mvf.susie.alpha:::.validate_multfsusie_response_controls(
      standardize,
      residual_variance,
      quantile_trans = TRUE
    ),
    "incompatible"
  )
})


test_that("mixed response standardization preserves modality scales", {
  Y <- list(
    Y_f = list(
      cbind(c(1, 3, 5), c(10, 14, 18)),
      cbind(c(2, 4, 6), c(3, 9, 15))
    ),
    Y_u = cbind(c(1, 3, 5), c(2, 6, 10))
  )

  observed <- multi_array_colScale(
    Y,
    scale = list(Y_f = c(FALSE, TRUE), Y_u = c(FALSE, TRUE))
  )

  expect_equal(observed$Y_u[, 1], c(-2, 0, 2))
  expect_equal(stats::sd(observed$Y_u[, 1]), 2)
  expect_equal(stats::sd(observed$Y_u[, 2]), 1)
  expect_equal(attr(observed$Y_u, "scaled:scale"), c(1, 4))

  expect_equal(colMeans(observed$Y_f[[1]]), c(0, 0))
  expect_equal(apply(observed$Y_f[[1]], 2, stats::sd), c(2, 4))
  expect_equal(mvf.susie.alpha:::.outcome_scale(observed$Y_f[[1]]), c(1, 1))
  expect_equal(apply(observed$Y_f[[2]], 2, stats::sd), c(1, 1))

  default_observed <- multi_array_colScale(Y)
  expect_equal(default_observed$Y_u, fsusieR::colScale(Y$Y_u))
  expect_equal(
    default_observed$Y_f,
    lapply(Y$Y_f, fsusieR::colScale)
  )
})


test_that("only requested residual variances are held fixed", {
  requested <- list(sd_f = c(0.4, NA), sd_u = c(NA, 2, NA))
  estimated <- list(sd_f = c(1.1, 1.2), sd_u = c(3.1, 3.2, 3.3))

  observed <- mvf.susie.alpha:::.apply_fixed_residual_variance(
    estimated,
    requested
  )
  expect_equal(observed$sd_f, c(0.4, 1.2))
  expect_equal(observed$sd_u, c(3.1, 2, 3.3))

  obj <- structure(
    list(n_wac = list(8L, 16L), n_univ = 3L),
    class = "multfsusie"
  )
  expect_equal(update_residual_variance(obj, observed)$sigma2, observed)
  expect_error(
    update_residual_variance(
      obj,
      list(sd_f = c(0.4, 1.2), sd_u = c(1, 2))
    ),
    "one positive residual variance"
  )
})


test_that("multfsusie can mix fixed and estimated univariate variances", {
  set.seed(101)
  n <- 40
  X <- matrix(rnorm(n * 5), nrow = n)
  Y_u <- cbind(
    3 + 0.5 * X[, 1] + rnorm(n, sd = 0.5),
    -2 + 2 * X[, 2] + rnorm(n, sd = 3)
  )

  fit <- multfsusie(
    Y = list(Y_f = NULL, Y_u = Y_u),
    X = X,
    L = 1,
    L_start = 1,
    verbose = FALSE,
    filter_cs = FALSE,
    posthoc = FALSE,
    standardize_Y = list(Y_u = c(FALSE, TRUE)),
    residual_variance = list(Y_u = c(0.25, NA))
  )

  expect_equal(fit$sigma2$sd_u[1], 0.25)
  expect_true(is.finite(fit$sigma2$sd_u[2]))
  expect_gt(fit$sigma2$sd_u[2], 0)
  expect_identical(fit$residual_variance_fixed$sd_u, c(TRUE, FALSE))
  expect_identical(fit$standardize_Y$Y_u, c(FALSE, TRUE))
  expect_equal(fit$csd_Y_u, c(1, stats::sd(Y_u[, 2])))
})


test_that("multfsusie can mix fixed and estimated functional variances", {
  set.seed(102)
  n <- 40
  X <- matrix(rnorm(n * 5), nrow = n)
  Y_f <- list(
    matrix(rnorm(n * 16, sd = 0.5), nrow = n) + X[, 1],
    matrix(rnorm(n * 16, sd = 2), nrow = n) + 0.5 * X[, 2]
  )

  fit <- multfsusie(
    Y = list(Y_f = Y_f, Y_u = NULL),
    X = X,
    L = 1,
    L_start = 1,
    verbose = FALSE,
    filter_cs = FALSE,
    posthoc = FALSE,
    post_processing = "none",
    standardize_Y = list(Y_f = c(FALSE, TRUE)),
    residual_variance = list(Y_f = c(0.5, NA))
  )

  expect_equal(fit$sigma2$sd_f[1], 0.5)
  expect_true(is.finite(fit$sigma2$sd_f[2]))
  expect_gt(fit$sigma2$sd_f[2], 0)
  expect_identical(fit$residual_variance_fixed$sd_f, c(TRUE, FALSE))
  expect_identical(fit$standardize_Y$Y_f, c(FALSE, TRUE))
  expect_equal(fit$csd_Y_f_pos[[1]], rep(1, 16))
  expect_equal(fit$csd_Y_f_wavelet[[1]], rep(1, 16))
  expect_length(fit$fitted_func[[1]], 2)
  expect_length(fit$fitted_func[[1]][[1]], 16)
})
