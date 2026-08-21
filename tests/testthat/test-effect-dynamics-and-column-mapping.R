make_dynamic_multfsusie <- function(L=3L, P=3L) {
  alpha <- lapply(seq_len(L), function(l) {
    value <- seq_len(P) + l
    value / sum(value)
  })
  structure(
    list(
      L=L,
      L_max=max(L, P),
      P=P,
      alpha=alpha,
      lBF=lapply(seq_len(L), function(l) seq_len(P) + l),
      lBF_per_trait=lapply(seq_len(L), function(l) {
        list(u_logBF=matrix(seq_len(P) + l, nrow=1))
      }),
      cs=lapply(seq_len(L), function(l) as.integer((l - 1L) %% P + 1L)),
      est_pi=lapply(seq_len(L), function(l) {
        list(est_pi_u=list(c(0.5, 0.5)), est_pi_f=NULL)
      }),
      est_sd=list(),
      fitted_u=lapply(seq_len(L), function(l) matrix(l, nrow=P, ncol=1)),
      fitted_u2=lapply(seq_len(L), function(l) matrix(l + 1, nrow=P, ncol=1)),
      fitted_wc=NULL,
      fitted_wc2=NULL,
      G_prior=list(
        G_prior_u=list(structure(list(), class="mixture_normal")),
        G_prior_f=NULL
      ),
      n_univ=1L,
      n_wac=NULL,
      tol_null_prior=0.001,
      pip=rep(0, P),
      KL=seq_len(L),
      lfsr_u=lapply(seq_len(L), function(l) l),
      lfsr_wc=NULL,
      alpha_hist=list(alpha),
      ind_fitted_val=list(),
      mean_X=seq_len(P),
      csd_X=rep(1, P),
      n_expand=0L,
      greedy=TRUE,
      greedy_backfit_update=FALSE,
      ELBO=numeric()
    ),
    class="multfsusie"
  )
}


test_that("multfsusie discard keeps every effect-indexed field aligned", {
  obj <- make_dynamic_multfsusie()
  obj$HMM_lBF <- lapply(seq_len(obj$L), function(l) list(l))
  retained_alpha <- obj$alpha[c(1, 3)]

  result <- discard_cs(obj, cs=2L, out_prep=TRUE)

  expect_equal(result$L, 2L)
  expect_equal(result$alpha, retained_alpha)
  for (field in c(
    "alpha", "lBF", "lBF_per_trait", "cs", "est_pi",
    "fitted_u", "fitted_u2", "KL", "lfsr_u", "HMM_lBF"
  )) {
    expect_length(result[[field]], result$L)
  }
  expect_equal(result$KL, c(1L, 3L))
  expect_length(result$pip, result$P)

  protected <- discard_cs(obj, cs=seq_len(obj$L), out_prep=TRUE)
  expect_equal(protected$L, 1L)
  expect_length(protected$cs, 1L)

  expect_error(discard_cs(obj, cs=0L), "invalid effect index")
  expect_error(discard_cs(obj, cs=4L), "invalid effect index")
})


test_that("multfsusie expansion respects P and keeps dynamic fields aligned", {
  obj <- make_dynamic_multfsusie(L=1L, P=3L)
  obj$L_max <- 5L

  expanded <- expand_multfsusie_obj(obj, L_extra=10L)

  expect_equal(expanded$L, 3L)
  for (field in c(
    "alpha", "lBF", "lBF_per_trait", "cs", "est_pi",
    "fitted_u", "fitted_u2", "KL", "lfsr_u"
  )) {
    expect_length(expanded[[field]], expanded$L)
  }
  expect_true(all(lengths(expanded$alpha) == expanded$P))
  expect_equal(expand_multfsusie_obj(expanded, 2L)$L, 3L)
  expect_error(expand_multfsusie_obj(obj, -1L), "non-negative integer")

  finalized <- obj
  finalized$column_index_space <- "original"
  finalized$fitted_P <- finalized$P
  expect_error(expand_multfsusie_obj(finalized, 1L), "after original-column")
})


test_that("multfsusie maps dynamic output without recalculating CS membership", {
  obj <- make_dynamic_multfsusie()
  obj$cs <- list(c(1L, 3L), c(1L, 2L), 2L)
  obj <- discard_cs(obj, cs=2L, out_prep=TRUE)
  obj$fitted_u <- lapply(seq_len(obj$L), function(l) as.numeric(l))
  obj$posthoc <- lapply(seq_len(obj$L), function(l) {
    list(variant_prior=c(0.2, 0.3, 0.5))
  })
  obj$purity <- as.list(rep(1, obj$L))

  restored <- mvf.susie.alpha:::restore_original_columns_multfsusie(
    obj,
    kept_index=c(1L, 3L, 5L),
    original_P=5L,
    removed_index=c(2L, 4L),
    variable_names=paste0("snp", 1:5),
    original_mean_X=10:14
  )

  expect_equal(restored$L, 2L)
  expect_equal(restored$P, 5L)
  expect_equal(restored$fitted_P, 3L)
  expect_equal(restored$variable_index, c(1L, 3L, 5L))
  expect_equal(restored$removed_variable_index, c(2L, 4L))
  expect_equal(
    restored$cs[[1]],
    structure(c(1L, 5L), names=c("snp1", "snp5"))
  )
  expect_equal(
    restored$cs[[2]],
    structure(3L, names="snp3")
  )
  expect_true(all(vapply(restored$alpha, length, integer(1)) == 5L))
  expect_true(all(vapply(restored$lBF, length, integer(1)) == 5L))
  expect_equal(restored$fitted_u, list(1, 2))
  expect_true(all(vapply(restored$fitted_u2, nrow, integer(1)) == 5L))
  expect_true(all(vapply(
    restored$lBF_per_trait,
    function(x) ncol(x$u_logBF),
    integer(1)
  ) == 5L))
  expect_equal(unname(restored$pip[c(2, 4)]), c(0, 0))
  expect_true(all(vapply(
    restored$posthoc,
    function(x) length(x$variant_prior),
    integer(1)
  ) == 5L))
  expect_equal(restored$mean_X, structure(10:14, names=paste0("snp", 1:5)))
  expect_equal(restored$n_cs, restored$L)
  expect_equal(restored$cs_size, c(2L, 1L))
  expect_identical(
    mvf.susie.alpha:::restore_original_columns_multfsusie(restored),
    restored
  )
})


test_that("postprocessed mixed-trait output handles constant and complete X", {
  make_postprocessed_output <- function(P) {
    obj <- make_dynamic_multfsusie(L=2L, P=P)
    obj$fitted_u <- list(c(1, 2, 3), c(4, 5, 6))
    obj$fitted_u2 <- lapply(seq_len(obj$L), function(l) {
      matrix(l, nrow=P, ncol=3L)
    })
    obj$fitted_wc <- lapply(seq_len(obj$L), function(l) {
      list(
        matrix(l, nrow=P, ncol=32L),
        matrix(l, nrow=P, ncol=64L)
      )
    })
    obj$fitted_wc2 <- lapply(seq_len(obj$L), function(l) {
      list(
        matrix(l + 1, nrow=P, ncol=32L),
        matrix(l + 1, nrow=P, ncol=64L)
      )
    })
    obj$lBF_per_trait <- lapply(seq_len(obj$L), function(l) {
      list(
        u_logBF=matrix(l, nrow=3L, ncol=P),
        f_logBF=matrix(l, nrow=2L, ncol=P)
      )
    })
    obj$cs <- list(c(3L, 4L), P)
    obj
  }

  with_constant <- make_postprocessed_output(P=99L)
  kept_index <- setdiff(seq_len(100L), 4L)
  restored <- NULL
  expect_error(
    restored <- mvf.susie.alpha:::restore_original_columns_multfsusie(
      with_constant,
      kept_index=kept_index,
      original_P=100L,
      removed_index=4L
    ),
    NA
  )

  # update_cal_fit_u() has already collapsed these to one value per
  # univariate trait; they must not be interpreted as SNP-indexed matrices.
  expect_identical(restored$fitted_u, with_constant$fitted_u)
  expect_true(all(vapply(restored$fitted_u2, nrow, integer(1)) == 100L))
  expect_true(all(vapply(restored$fitted_wc, function(effect) {
    all(vapply(effect, nrow, integer(1)) == 100L)
  }, logical(1))))
  expect_true(all(vapply(restored$fitted_wc2, function(effect) {
    all(vapply(effect, nrow, integer(1)) == 100L)
  }, logical(1))))
  expect_equal(restored$cs[[1]], c(3L, 5L))
  expect_equal(restored$cs[[2]], 100L)
  expect_equal(restored$pip[4L], 0)

  without_constant <- make_postprocessed_output(P=100L)
  unchanged_columns <- NULL
  expect_error(
    unchanged_columns <- mvf.susie.alpha:::restore_original_columns_multfsusie(
      without_constant,
      kept_index=seq_len(100L),
      original_P=100L
    ),
    NA
  )
  expect_identical(unchanged_columns$fitted_u, without_constant$fitted_u)
  expect_equal(unchanged_columns$cs, without_constant$cs)
})


test_that("posthoc NULL effects retain their effect positions", {
  obj <- make_dynamic_multfsusie(L=2L, P=3L)
  obj$alpha[[1]] <- rep(0, 3L)

  result <- mvf.susie.alpha:::posthoc_multfsusie(obj)

  expect_length(result, 2L)
  expect_null(result[[1]])
  expect_true(is.list(result[[2]]))
  expect_true("posthoc" %in% names(result[[2]]))
})


test_that("constant-column detection covers modality subsets and one-column X", {
  X <- cbind(
    c(0, 1, 2, 3),
    rep(1, 4),
    c(0, 0, 1, 2),
    c(0, 1, 0, 1)
  )
  groups <- list(idx_f=list(1:2), idx_u=list(1:4))
  expect_equal(
    mvf.susie.alpha:::check_cst_X_sub_case(X, groups),
    c(2L, 3L)
  )
  expect_length(
    mvf.susie.alpha:::check_cst_X_sub_case(
      matrix(1:4, ncol=1),
      list(idx_f=list(1:4), idx_u=NULL)
    ),
    0L
  )
  expect_equal(
    mvf.susie.alpha:::check_cst_X_sub_case(
      matrix(1, nrow=4, ncol=1),
      list(idx_f=list(1:4), idx_u=NULL)
    ),
    1L
  )
  expect_error(
    mvf.susie.alpha:::check_cst_X_sub_case(
      matrix(1:4, ncol=1),
      list(idx_f=list(1L), idx_u=NULL)
    ),
    "At least two"
  )
})

test_that("singleton credible sets are not rejected by the purity criterion", {
  obj <- make_dynamic_multfsusie(L=2L, P=2L)
  obj$cs <- list(1L, 2L)
  obj$est_pi <- lapply(seq_len(obj$L), function(l) {
    list(est_pi_u=list(c(0.5, 0.5)), est_pi_f=NULL)
  })
  obj$tol_null_prior <- 0.001

  expect_length(
    which_dummy_cs(obj, X=matrix(c(0, 1, 0, 1), nrow=2)),
    0L
  )
})
