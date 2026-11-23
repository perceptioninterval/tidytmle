test_that("Automatic informative missingness modeling triggers when outcome has NAs", {
  df <- sim_causal_cohort(n = 350, outcome_type = "continuous", censoring_rate = 0.2, seed = 606)
  # Drop explicit delta column so tidytmle auto-detects NAs in y
  df_no_delta <- dplyr::select(df, -delta)

  expect_message(
    fit <- tmle_spec(
      outcome = y,
      treatment = a,
      covariates = c(age, sex, starts_with("lab_"))
    ) |>
      set_resamples(v_folds = 4L) |>
      fit(data = df_no_delta),
    "missingness indicator"
  )

  expect_true(fit$cache$has_censoring)
  expect_true(fit$cache$auto_censoring)
  expect_false(is.null(fit$nuisance$C))
  gl <- glance(fit)
  expect_true(gl$n_observed < gl$n)
  expect_true(abs(gl$mean_eif) < 1e-4)
})

test_that("Explicit censoring indicator works cleanly", {
  df <- sim_causal_cohort(n = 350, outcome_type = "continuous", censoring_rate = 0.15, seed = 607)

  fit <- tmle_spec(
    outcome = y,
    treatment = a,
    covariates = c(age, sex, starts_with("lab_")),
    censoring = delta
  ) |>
    set_resamples(v_folds = 4L) |>
    fit(data = df)

  expect_true(fit$cache$has_censoring)
  expect_false(fit$cache$auto_censoring)
  expect_true("C" %in% extract_sl_weights(fit)$nuisance)
  expect_true(abs(glance(fit)$mean_eif) < 1e-4)
})

test_that("Clustered data triggers grouped CV and cluster-robust EIF standard errors", {
  df <- sim_causal_cohort(n = 400, outcome_type = "continuous", n_clusters = 20L, seed = 707)

  fit_cl <- tmle_spec(
    outcome = y,
    treatment = a,
    covariates = c(age, sex, starts_with("lab_")),
    id = clinic_id
  ) |>
    set_resamples(v_folds = 5L) |>
    fit(data = df)

  gl <- glance(fit_cl)
  expect_equal(gl$n_clusters, 20L)
  expect_true(abs(gl$mean_eif) < 1e-4)
  td <- tidy(fit_cl)
  expect_true(is.finite(td$std.error) && td$std.error > 0)
})

test_that("Native hardhat importance_weights and frequency_weights are supported", {
  df <- sim_causal_cohort(n = 300, outcome_type = "continuous", seed = 808)
  df$imp_wt <- hardhat::importance_weights(runif(nrow(df), 0.5, 2.0))
  df$freq_wt <- hardhat::frequency_weights(sample(1:3, nrow(df), replace = TRUE))

  fit_imp <- tmle_spec(
    outcome = y,
    treatment = a,
    covariates = c(age, sex, starts_with("lab_")),
    weights = imp_wt
  ) |>
    set_resamples(v_folds = 4L) |>
    fit(data = df)

  fit_freq <- tmle_spec(
    outcome = y,
    treatment = a,
    covariates = c(age, sex, starts_with("lab_")),
    weights = freq_wt
  ) |>
    set_resamples(v_folds = 4L) |>
    fit(data = df)

  expect_s3_class(fit_imp, "tmle_fit")
  expect_s3_class(fit_freq, "tmle_fit")
  expect_true(abs(glance(fit_imp)$mean_eif) < 1e-4)
  expect_true(abs(glance(fit_freq)$mean_eif) < 1e-4)

  # Extreme weights must genuinely alter both nuisance predictions and estimates
  df$extreme_wt <- hardhat::importance_weights(ifelse(df$age > 60, 10.0, 0.2))
  fit_unw <- tmle_spec(y, a, c(age, sex, lab_1, lab_2, lab_3)) |>
    set_resamples(v_folds = 3L) |>
    fit(data = df)
  fit_ext <- tmle_spec(y, a, c(age, sex, lab_1, lab_2, lab_3), weights = extreme_wt) |>
    set_resamples(v_folds = 3L) |>
    fit(data = df)
  expect_true(abs(tidy(fit_unw)$estimate - tidy(fit_ext)$estimate) > 1e-4)
})

test_that("Missing baseline covariates are auto-imputed cleanly at Level 1", {
  df <- sim_causal_cohort(n = 250, outcome_type = "continuous", seed = 888)
  df$lab_1[1:15] <- NA_real_
  df$sex[16:25] <- NA

  fit <- tmle_spec(y, a, c(age, sex, starts_with("lab_"))) |>
    set_resamples(v_folds = 3L) |>
    fit(data = df)

  expect_s3_class(fit, "tmle_fit")
  expect_equal(glance(fit)$n, 250L)
  expect_true(abs(glance(fit)$mean_eif) < 1e-4)
})

test_that("Ratio guardrails catch invalid continuous OR / non-positive RR, and ATT shift is exact", {
  df_cont <- sim_causal_cohort(n = 150, outcome_type = "continuous", seed = 919)
  expect_error(
    tmle_spec(y, a, c(age, sex, lab_1), estimand = "or") |>
      set_resamples(v_folds = 3L) |>
      fit(data = df_cont),
    "Odds Ratio"
  )

  df_neg <- df_cont
  df_neg$y <- df_neg$y - 10
  expect_error(
    tmle_spec(y, a, c(age, sex, lab_1), estimand = "rr") |>
      set_resamples(v_folds = 3L) |>
      fit(data = df_neg),
    "strictly positive"
  )

  fit_att <- tmle_spec(y, a, c(age, sex, lab_1, lab_2), estimand = "att") |>
    set_resamples(v_folds = 3L) |>
    fit(data = df_cont)
  gl <- glance(fit_att)
  expect_true(abs(gl$mean_eif) < 1e-4)
  plug_att <- mean((fit_att$targeting$q1_init - fit_att$targeting$q0_init)[df_cont$a == 1L])
  expect_equal(gl$plugin_vs_tmle_shift, tidy(fit_att)$estimate - plug_att, tolerance = 1e-6)
})

test_that("Bounded logistic fluctuation keeps continuous predictions in observed range under positivity stress", {
  df <- sim_causal_cohort(n = 400, outcome_type = "continuous", positivity_stress = TRUE, seed = 909)

  fit <- tmle_spec(
    outcome = y,
    treatment = a,
    covariates = c(age, sex, starts_with("lab_"))
  ) |>
    set_resamples(v_folds = 4L) |>
    fit(data = df)

  aug <- augment(fit)
  y_rng <- range(df$y, na.rm = TRUE)
  expect_true(all(aug$.q0_star >= y_rng[1] - 1e-8 & aug$.q0_star <= y_rng[2] + 1e-8))
  expect_true(all(aug$.q1_star >= y_rng[1] - 1e-8 & aug$.q1_star <= y_rng[2] + 1e-8))
  expect_true(glance(fit)$pct_truncated > 0)
  expect_true(abs(glance(fit)$mean_eif) < 1e-4)
})
