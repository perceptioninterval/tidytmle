test_that("Modified Treatment Policies shift() and scale_tx() work for continuous treatments", {
  df_shift <- sim_causal_cohort(n = 400, outcome_type = "continuous", treatment_type = "continuous", seed = 1212)
  truth <- attr(df_shift, "truth")

  fit_s <- tmle_spec(
    outcome = y,
    treatment = a,
    covariates = c(age, sex, starts_with("lab_")),
    estimand = shift(1.0)
  ) |>
    set_resamples(v_folds = 4L) |>
    fit(data = df_shift)

  td_s <- tidy(fit_s)
  expect_equal(td_s$estimand, "shift(1)")
  expect_true(abs(td_s$estimate - truth$true_shift_1) < 3.5 * td_s$std.error)
  expect_true(abs(glance(fit_s)$mean_eif) < 1e-4)

  fit_sc <- tmle_spec(
    outcome = y,
    treatment = a,
    covariates = c(age, sex, starts_with("lab_")),
    estimand = scale_tx(0.9)
  ) |>
    set_resamples(v_folds = 4L) |>
    fit(data = df_shift)

  td_sc <- tidy(fit_sc)
  expect_equal(td_sc$estimand, "scale_tx(0.9)")
  expect_true(abs(td_sc$estimate - truth$true_scale_09) < 3.5 * td_sc$std.error)
  expect_true(abs(glance(fit_sc)$mean_eif) < 1e-4)
})

test_that("DR-Learner CATE estimation (fit_cate, predict, tidy, augment, qini_auc) recovers heterogeneity", {
  df <- sim_causal_cohort(n = 500, outcome_type = "continuous", seed = 1313)
  truth <- attr(df, "truth")

  fit <- tmle_spec(
    outcome = y,
    treatment = a,
    covariates = c(age, sex, starts_with("lab_")),
    estimand = "ate"
  ) |>
    set_resamples(v_folds = 4L) |>
    fit(data = df)

  cate_fit <- fit_cate(fit, covariates = c(age, sex, lab_1, lab_2, lab_3), v_folds = 4L)
  expect_s3_class(cate_fit, "tmle_cate")
  expect_output(print(cate_fit), "Doubly Robust CATE")

  preds <- predict(cate_fit, new_data = df[1:25, ])
  expect_equal(nrow(preds), 25L)
  expect_true(".pred" %in% names(preds))

  aug_c <- augment(cate_fit)
  expect_true(all(c(".cate", ".eif_pseudo") %in% names(aug_c)))
  # Estimated CATE should positively correlate with true individual treatment effects
  expect_true(stats::cor(aug_c$.cate, truth$true_cate) > 0.35)

  blp <- tidy(cate_fit)
  expect_equal(nrow(blp), 2L)
  q_tbl <- qini_auc(cate_fit)
  expect_equal(q_tbl$.metric, c("qini_auc", "auuc"))
  expect_true(q_tbl$.estimate[1L] > 0)

  expect_s3_class(autoplot(cate_fit, type = "qini"), "ggplot")
  expect_s3_class(autoplot(cate_fit, type = "cate"), "ggplot")
})
