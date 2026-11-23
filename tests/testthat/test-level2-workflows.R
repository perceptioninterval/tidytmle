test_that("Level 2 supports single parsnip spec, workflow, and workflow_set learners", {
  df <- sim_causal_cohort(n = 300, outcome_type = "continuous", seed = 404)

  base_rec <- recipes::recipe(y ~ ., data = df) |>
    recipes::step_impute_median(recipes::all_numeric_predictors()) |>
    recipes::step_dummy(recipes::all_nominal_predictors()) |>
    recipes::step_zv(recipes::all_predictors())

  has_glmnet <- requireNamespace("glmnet", quietly = TRUE)
  q_models <- list(
    lm = parsnip::linear_reg() |> parsnip::set_engine("glm")
  )
  if (has_glmnet) {
    q_models$enet <- parsnip::linear_reg(penalty = 0.05, mixture = 0.5) |>
      parsnip::set_engine("glmnet")
  } else {
    q_models$lm_alt <- parsnip::linear_reg() |> parsnip::set_engine("lm")
  }

  q_wf_set <- workflowsets::workflow_set(
    preproc = list(base = base_rec),
    models = q_models
  )

  g_single <- parsnip::logistic_reg() |> parsnip::set_engine("glm")

  fit <- tmle_spec(
    outcome = y,
    treatment = a,
    covariates = c(age, sex, starts_with("lab_")),
    estimand = "ate"
  ) |>
    set_preprocessing(base_rec) |>
    set_outcome_learners(q_wf_set) |>
    set_treatment_learners(g_single) |>
    set_resamples(v_folds = 4L) |>
    fit(data = df)

  expect_s3_class(fit, "tmle_fit")
  w_tbl <- extract_sl_weights(fit)
  expect_true("base_lm" %in% w_tbl$learner[w_tbl$nuisance == "Q"])
  if (has_glmnet) {
    expect_true("base_enet" %in% w_tbl$learner[w_tbl$nuisance == "Q"])
  } else {
    expect_true("base_lm_alt" %in% w_tbl$learner[w_tbl$nuisance == "Q"])
  }
  expect_equal(sum(w_tbl$weight[w_tbl$nuisance == "Q"]), 1, tolerance = 1e-6)
  expect_equal(sum(w_tbl$weight[w_tbl$nuisance == "g"]), 1, tolerance = 1e-6)
})

test_that("Role-aware recipe splitting strictly prevents treatment and outcome leakage into g", {
  df <- sim_causal_cohort(n = 200, outcome_type = "continuous", seed = 505)
  df$clinic_id <- factor(sample(paste0("c", 1:10), nrow(df), replace = TRUE))

  # User writes a naive recipe y ~ . containing treatment 'a' and cluster 'clinic_id'
  base_rec <- recipes::recipe(y ~ ., data = df) |>
    recipes::step_dummy(recipes::all_nominal_predictors()) |>
    recipes::step_zv(recipes::all_predictors())

  spec <- tmle_spec(
    outcome = y,
    treatment = a,
    covariates = c(age, sex, lab_1, lab_2, lab_3),
    id = clinic_id
  ) |>
    set_preprocessing(base_rec)

  roles <- tidytmle:::resolve_causal_roles(spec, df)
  ndata <- tidytmle:::build_nuisance_data(roles)
  recs <- tidytmle:::split_causal_recipes(base_rec, roles, ndata)

  # Bake prepped recipes and verify exact columns
  prepped_Q <- recipes::prep(recs$Q, training = ndata$Q)
  baked_Q <- recipes::bake(prepped_Q, new_data = NULL)
  expect_true("a" %in% names(baked_Q))
  expect_true("y" %in% names(baked_Q))
  expect_false(any(grepl("^clinic_id", names(baked_Q))))

  prepped_g <- recipes::prep(recs$g, training = ndata$g)
  baked_g <- recipes::bake(prepped_g, new_data = NULL)
  # 'a' must be outcome, NEVER predictor, and 'y' and 'clinic_id' must be excluded!
  g_preds <- setdiff(names(baked_g), "a")
  expect_false("y" %in% g_preds)
  expect_false("a" %in% g_preds)
  expect_false(any(grepl("^clinic_id", g_preds)))
  expect_true(all(c("age", "lab_1", "lab_2", "lab_3") %in% g_preds))
})
