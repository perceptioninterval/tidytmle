test_that("Level 1 3-line safe default works on continuous outcomes and recovers truth", {
  df <- sim_causal_cohort(n = 400, outcome_type = "continuous", seed = 101)
  truth <- attr(df, "truth")

  fit <- tmle_spec(
    outcome = y,
    treatment = a,
    covariates = c(age, sex, starts_with("lab_")),
    estimand = "ate"
  ) |>
    set_resamples(v_folds = 5L) |>
    fit(data = df)

  expect_s3_class(fit, "tmle_fit")
  td <- tidy(fit)
  expect_s3_class(td, "tbl_df")
  expect_equal(td$estimand, "ate")
  expect_true(is.finite(td$estimate) && is.finite(td$std.error) && td$std.error > 0)
  # Unbiased 95% CI should cover true ATE (within reasonable sampling variance)
  expect_true(abs(td$estimate - truth$true_ate) < 3.5 * td$std.error)

  gl <- glance(fit)
  expect_equal(nrow(gl), 1L)
  expect_true(abs(gl$mean_eif) < 1e-4)
  expect_true(gl$cv_tmle)
  expect_equal(gl$v_folds, 5L)
})

test_that("Level 1 supports simultaneous binary multi-estimand output (ATE, RR, OR, EY1, EY0)", {
  df <- sim_causal_cohort(n = 450, outcome_type = "binary", seed = 202)
  truth <- attr(df, "truth")

  fit <- tmle_spec(
    outcome = y,
    treatment = a,
    covariates = c(age, sex, starts_with("lab_")),
    estimand = c("ate", "rr", "or", "ey1", "ey0")
  ) |>
    set_resamples(v_folds = 5L) |>
    fit(data = df)

  td <- tidy(fit)
  expect_equal(td$estimand, c("ate", "rr", "or", "ey1", "ey0"))
  # RR and OR CIs must stay strictly positive (log scale Delta method)
  rr_row <- td[td$estimand == "rr", ]
  or_row <- td[td$estimand == "or", ]
  expect_true(rr_row$conf.low > 0 && rr_row$conf.high > rr_row$estimate)
  expect_true(or_row$conf.low > 0 && or_row$conf.high > or_row$estimate)

  # Joint targeting guarantees near-zero EIF for both arms
  expect_true(abs(mean(fit$targeting$eif_1)) < 1e-4)
  expect_true(abs(mean(fit$targeting$eif_0)) < 1e-4)
  expect_true(abs(td$estimate[td$estimand == "ate"] - truth$true_ate) < 0.20)
})

test_that("Level 1 works with omitted covariates argument and factor/logical treatment & outcome", {
  df <- sim_causal_cohort(n = 250, outcome_type = "binary", seed = 303)
  df$a <- factor(ifelse(df$a == 1, "treated", "control"), levels = c("control", "treated"))
  df$y <- as.logical(df$y)

  fit <- tmle_spec(outcome = y, treatment = a) |>
    set_resamples(v_folds = 4L) |>
    fit(data = df)

  expect_s3_class(fit, "tmle_fit")
  expect_equal(sort(fit$cache$covariate_cols), sort(c("age", "sex", "lab_1", "lab_2", "lab_3")))
  expect_true(abs(glance(fit)$mean_eif) < 1e-4)
})
