test_that("Instant retarget() updates g_bounds, estimands (ATT/ATC/RR/OR), and subgroups without refitting", {
  df <- sim_causal_cohort(n = 400, outcome_type = "binary", positivity_stress = TRUE, seed = 1010)

  fit <- tmle_spec(
    outcome = y,
    treatment = a,
    covariates = c(age, sex, starts_with("lab_")),
    estimand = "ate"
  ) |>
    set_resamples(v_folds = 4L) |>
    fit(data = df)

  t0 <- proc.time()[["elapsed"]]
  fit_strict <- retarget(fit, g_bounds = c(0.10, 0.90))
  fit_ratios <- retarget(fit, estimand = c("ate", "rr", "or", "att", "atc"))
  fit_older <- retarget(fit, subgroup = age >= 55)
  dt <- proc.time()[["elapsed"]] - t0

  expect_true(dt < 1.0)
  expect_equal(glance(fit_strict)$g_bound_low, 0.10)
  expect_equal(glance(fit_strict)$g_bound_high, 0.90)
  expect_true(glance(fit_strict)$pct_truncated >= glance(fit)$pct_truncated)

  td_r <- tidy(fit_ratios)
  expect_equal(td_r$estimand, c("ate", "rr", "or", "att", "atc"))
  expect_true(all(abs(fit_ratios$targeting$mean_eif) < 1e-4))

  gl_old <- glance(fit_older)
  expect_equal(gl_old$n, sum(df$age >= 55))
  expect_true(abs(gl_old$mean_eif) < 1e-4)
})

test_that("All S3 verbs (tidy simultaneous/e_value, glance, augment, collect_metrics, extract_sl_weights, print, autoplot) work", {
  df <- sim_causal_cohort(n = 300, outcome_type = "binary", seed = 1111)

  fit <- tmle_spec(
    outcome = y,
    treatment = a,
    covariates = c(age, sex, starts_with("lab_")),
    estimand = c("ate", "rr", "or")
  ) |>
    set_resamples(v_folds = 4L) |>
    fit(data = df)

  # tidy with simultaneous uniform bands and E-values
  td_full <- tidy(fit, simultaneous = TRUE, e_value = TRUE)
  expect_true(all(c("e_value", "e_value_ci", "crit.value") %in% names(td_full)))
  expect_true(all(td_full$crit.value >= stats::qnorm(0.975) - 0.05))
  expect_true(all(td_full$e_value >= 1))

  # glance
  gl <- glance(fit)
  expect_equal(nrow(gl), 1L)
  expect_true(all(c("ess_treated", "ess_control", "max_eif_leverage") %in% names(gl)))

  # augment
  aug <- augment(fit)
  expect_equal(nrow(aug), nrow(df))
  expect_true(all(c(
    ".q0_init", ".q1_init", ".q0_star", ".q1_star",
    ".g_hat", ".g_trunc", ".clever_h", ".ite_star", ".eif", ".eif_leverage"
  ) %in% names(aug)))
  expect_equal(sum(aug$.eif_leverage), 1, tolerance = 1e-6)

  # collect_metrics
  mets_sum <- collect_metrics(fit, summarize = TRUE)
  mets_fold <- collect_metrics(fit, summarize = FALSE)
  expect_true("SuperLearner" %in% mets_sum$learner)
  expect_true(all(c("Q", "g") %in% mets_sum$nuisance))
  expect_true("fold" %in% names(mets_fold))

  # extract_sl_weights
  sl_w <- extract_sl_weights(fit)
  sl_w_fold <- extract_sl_weights(fit, by_fold = TRUE)
  expect_true(all(c("nuisance", "learner", "weight") %in% names(sl_w)))
  expect_true("fold" %in% names(sl_w_fold))

  # print methods
  expect_output(print(fit$spec), "Targeted Maximum Likelihood")
  expect_output(print(fit), "TMLE Fit Summary")

  # All 6 autoplot diagnostics
  for (pt in c("overlap", "targeting", "learners", "eif", "sensitivity", "balance")) {
    p <- autoplot(fit, type = pt, g_grid = c(0.02, 0.05, 0.10))
    expect_s3_class(p, "ggplot")
  }
})

test_that("save_nuisance policies and mirai fine-grained parallelism work cleanly", {
  df <- sim_causal_cohort(n = 180, outcome_type = "continuous", seed = 1414)

  fit_pred <- tmle_spec(y, a, c(age, sex, lab_1, lab_2)) |>
    set_resamples(v_folds = 3L) |>
    set_targeting(save_nuisance = "predictions") |>
    fit(data = df)
  expect_null(fit_pred$nuisance$Q$models)

  fit_butcher <- tmle_spec(y, a, c(age, sex, lab_1, lab_2)) |>
    set_resamples(v_folds = 3L) |>
    set_targeting(save_nuisance = "butcher") |>
    fit(data = df)
  expect_true(length(fit_butcher$nuisance$Q$models) > 0L)

  if (requireNamespace("mirai", quietly = TRUE)) {
    mirai::daemons(2L)
    on.exit(mirai::daemons(0L), add = TRUE)
    fit_par <- tmle_spec(y, a, c(age, sex, lab_1, lab_2)) |>
      set_resamples(v_folds = 3L) |>
      set_targeting(parallel = TRUE) |>
      fit(data = df)
    expect_s3_class(fit_par, "tmle_fit")
    expect_true(abs(glance(fit_par)$mean_eif) < 1e-4)
  }
})
