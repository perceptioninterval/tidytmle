# tidytmle

**Cross-Fitted Targeted Maximum Likelihood Estimation for `tidymodels`**

`tidytmle` is a semiparametric causal inference orchestration layer built natively on top of `workflows`, `workflowsets`, `recipes`, `rsample`, `hardhat`, and `yardstick`. Instead of forcing TMLE into a prediction-only `parsnip` model specification, `tidytmle` composes any `tidymodels` workflow or `workflow_set` as cross-fitted Super Learner nuisance models while enforcing statistical guardrails automatically.

## Key Capabilities

- **Progressive-Disclosure Specification API**: Start with a 3-line Level 1 specification (`tmle_spec(outcome, treatment, covariates, estimand) |> fit(data)`) and scale seamlessly to Level 2 `recipes`, `parsnip` specs, `workflows`, and `workflow_set` nuisance ensembles.
- **Role-Aware Recipe Cloning**: Automatically splits a base `recipe` across $Q(A, W)$, $g(W)$, and $C(A, W)$ so treatment $A$ never leaks into propensity models and `id`/censoring/weight columns never leak into predictors.
- **Automatic Informative-Missingness Handling**: Detects `NA`s in `outcome` when `censoring` is omitted and introduces $\Delta = \mathbb{I}(\neg\text{is.na}(Y))$ modeled as $C(A, W) = P(\Delta = 1 \mid A, W)$.
- **Built-In Semiparametric Guardrails**: 10-fold CV-TMLE enabled by default, bounded canonical logistic fluctuation on scaled `[0, 1]` continuous outcomes, adaptive sample-size positivity bounds (`5 / (sqrt(n) * log(n))`), joint EIF targeting for simultaneous ATE / RR / OR inference with log-scale CIs, cluster-robust EIF standard errors, native `hardhat::importance_weights()` and `frequency_weights()`, simultaneous multiplier bootstrap bands, and E-values.
- **Instant `retarget()`**: Re-run positivity bounds, estimands, or subgroup analyses in milliseconds using cached out-of-fold nuisance predictions with zero ML refitting.
- **Modified Treatment Policies & DR-Learner CATE**: Built-in `shift()` / `scale_tx()` density-ratio TMLE for continuous exposures and `fit_cate()` + `qini_auc()` for heterogeneous treatment effect estimation.
- **Standard Tidymodels Verbs**: `tidy()`, `glance()`, `augment()`, `collect_metrics()`, `extract_sl_weights()`, `autoplot()` (`"overlap"`, `"learners"`, `"eif"`, `"sensitivity"`, `"balance"`), and `cli` diagnostic summaries.

## Quick Example

```r
library(tidytmle)

df <- sim_causal_cohort(n = 500, outcome_type = "binary", seed = 42)

# Level 1: 3-line safe default
fit <- tmle_spec(
  outcome    = y,
  treatment  = a,
  covariates = c(age, sex, starts_with("lab_")),
  estimand   = c("ate", "rr", "or")
) |>
  fit(data = df)

tidy(fit, e_value = TRUE)
glance(fit)

# Instant retarget() with zero ML refitting
fit_strict <- retarget(fit, g_bounds = c(0.05, 0.95))
fit_older  <- retarget(fit, subgroup = age >= 60)
```
