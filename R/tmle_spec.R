#' Specify a Targeted Maximum Likelihood Estimation (TMLE) Pipeline
#'
#' Creates a progressive-disclosure specification for cross-fitted Targeted
#' Maximum Likelihood Estimation (CV-TMLE). At Level 1, only `outcome`,
#' `treatment`, `covariates`, and `estimand` are needed; sensible defaults for
#' 10-fold cross-fitting, Super Learning, adaptive positivity truncation, and
#' bounded logistic fluctuation are enabled automatically. At Level 2, attach
#' custom `recipes`, `parsnip` specs, `workflows`, `workflowsets`, and `rsample`
#' objects using the `set_*()` verbs.
#'
#' @details
#' ## How TMLE Works ("Predict, Then Target")
#' Standard machine learning regression models \eqn{\bar{Q}^0(A, W) = E[Y \mid A, W]}
#' optimize global prediction accuracy rather than the bias-variance trade-off
#' for a causal contrast such as the Average Treatment Effect (ATE). Conversely,
#' inverse-probability weighting (IPW) relies entirely on the propensity score
#' \eqn{g(W) = P(A = 1 \mid W)} and blows up in variance when positivity is
#' stressed.
#'
#' As synthesized in Schuler & Rose (2017) and Luque-Fernandez et al. (2018),
#' TMLE resolves both weaknesses in two stages:
#' 1. **Predict**: Estimate the initial outcome surface \eqn{\bar{Q}^0(A, W)}
#'    and treatment propensity \eqn{g(W)} (plus missingness/censoring
#'    \eqn{C(A, W)} when outcomes are missing) using cross-validated Super
#'    Learning.
#' 2. **Target (Fluctuate)**: Construct arm-specific Clever Covariates
#'    \eqn{H_1(A, W) = \mathbb{I}(A = 1) / g(W)} and
#'    \eqn{H_0(A, W) = \mathbb{I}(A = 0) / (1 - g(W))} (including censoring
#'    weights when outcomes are missing) derived from the Efficient Influence
#'    Function (EIF). By default (`fluctuation = "logistic"`), outcomes are
#'    mapped to `(0, 1)` and fluctuated via arm-specific weighted-intercept
#'    canonical logistic models \eqn{\text{logit}\,\bar{Q}^*(a, W) =
#'    \text{logit}\,\bar{Q}^0(a, W) + \epsilon_a} using weights \eqn{H_a}
#'    (matching `target.gwt = TRUE`). Classical unweighted Clever-Covariate
#'    fluctuation (\eqn{\text{logit}\,\bar{Q}^0 + \epsilon_1 H_1 + \epsilon_0 H_0},
#'    matching `target.gwt = FALSE`) is available via `fluctuation =
#'    "logistic_covariate"`. Both nudge initial predictions until the sample mean
#'    of each component EIF solves the zero score equation.
#'
#' This makes TMLE a **doubly robust**, **locally efficient substitution
#' estimator**: point estimation remains consistent if *either* \eqn{\bar{Q}}
#' *or* \eqn{g} (and \eqn{C} under missingness) is consistently estimated, while
#' EIF-based Wald confidence intervals and semiparametric efficiency hold when
#' *both* nuisance models converge fast enough that their second-order product
#' remainder is \eqn{o_p(n^{-1/2})}.
#'
#' @references
#' van der Laan, M. J., & Rubin, D. (2006). Targeted maximum likelihood learning.
#' *The International Journal of Biostatistics*, 2(1), Article 11.
#'
#' Schuler, M. S., & Rose, S. (2017). Targeted maximum likelihood estimation for
#' causal inference in observational studies. *American Journal of Epidemiology*,
#' 185(1), 65–73.
#'
#' Luque-Fernandez, M. A., Schomaker, M., Rachet, B., & Schnitzer, M. E.
#' (2018). Targeted maximum likelihood estimation for a binary treatment: A
#' tutorial. *Statistics in Medicine*, 37(16), 2530–2546.
#'
#' @param outcome Tidyselect specification of the outcome column. May be numeric
#'   (continuous or 0/1 binary) or a 2-level factor/logical.
#' @param treatment Tidyselect specification of the treatment/exposure column.
#'   For standard causal estimands (`"ate"`, `"att"`, `"atc"`, `"rr"`, `"or"`,
#'   `"ey1"`, `"ey0"`), must be binary (0/1, logical, or 2-level factor). For
#'   Modified Treatment Policies (`shift()`, `scale_tx()`), numeric continuous.
#' @param covariates Optional tidyselect specification of baseline confounders
#'   \eqn{W}. Defaults to `NULL`, which automatically selects all columns in
#'   `data` other than `outcome`, `treatment`, `censoring`, `id`, and `weights`.
#' @param censoring Optional tidyselect specification of the observation/non-censoring
#'   indicator \eqn{\Delta \in \{0, 1\}} (where 1 = observed, 0 = censored).
#'   If omitted and `outcome` contains `NA`s, `tidytmle` automatically
#'   introduces \eqn{\Delta = \mathbb{I}(\neg\text{is.na}(Y))} and models
#'   informative missingness via \eqn{C(A, W) = P(\Delta = 1 \mid A, W)}.
#' @param id Optional tidyselect specification of a clustering/grouping identifier.
#'   Automatically triggers grouped cross-validation and cluster-robust
#'   influence-function standard errors.
#' @param weights Optional tidyselect specification of a case-weights column
#'   (supporting `hardhat::importance_weights()`, `hardhat::frequency_weights()`,
#'   or positive numeric vectors).
#' @param estimand Causal estimand(s) to target: one or more of `"ate"`, `"att"`,
#'   `"atc"`, `"rr"`, `"or"`, `"ey1"`, `"ey0"`, or a Modified Treatment Policy
#'   created by [shift()] or [scale_tx()]. Defaults to `"ate"`.
#'
#' @return An object of class `"tmle_spec"`.
#'
#' @examples
#' df <- sim_causal_cohort(n = 150, outcome_type = "binary", seed = 42)
#' spec <- tmle_spec(
#'   outcome = y,
#'   treatment = a,
#'   covariates = c(age, sex, starts_with("lab_")),
#'   estimand = c("ate", "rr", "or")
#' )
#' spec
#'
#' @export
tmle_spec <- function(outcome,
                      treatment,
                      covariates = NULL,
                      censoring = NULL,
                      id = NULL,
                      weights = NULL,
                      estimand = "ate") {
  out_quo <- rlang::enquo(outcome)
  trt_quo <- rlang::enquo(treatment)
  cov_quo <- rlang::enquo(covariates)
  cen_quo <- rlang::enquo(censoring)
  id_quo <- rlang::enquo(id)
  wt_quo <- rlang::enquo(weights)

  if (rlang::quo_is_missing(out_quo) || rlang::quo_is_null(out_quo)) {
    cli::cli_abort("{.arg outcome} must be supplied.")
  }
  if (rlang::quo_is_missing(trt_quo) || rlang::quo_is_null(trt_quo)) {
    cli::cli_abort("{.arg treatment} must be supplied.")
  }

  est_norm <- normalize_estimand(estimand)

  structure(
    list(
      roles = list(
        outcome = out_quo,
        treatment = trt_quo,
        covariates = cov_quo,
        censoring = cen_quo,
        id = id_quo,
        weights = wt_quo
      ),
      estimand = est_norm,
      raw_estimand = estimand,
      preprocessing = NULL,
      learners = list(
        Q = list(spec = "default", metalearner = "auto", stratify = FALSE),
        g = list(spec = "default", metalearner = "auto"),
        C = list(spec = "default", metalearner = "auto")
      ),
      resamples = list(
        obj = NULL,
        v_folds = 10L,
        strata = TRUE
      ),
      targeting = list(
        cv_tmle = TRUE,
        g_bounds = "adaptive",
        fluctuation = "logistic",
        max_iter = 100L,
        tol = 1e-6,
        save_nuisance = "predictions",
        parallel = FALSE
      )
    ),
    class = "tmle_spec"
  )
}

check_is_tmle_spec <- function(spec, call = rlang::caller_env()) {
  if (!inherits(spec, "tmle_spec")) {
    cli::cli_abort(
      "Expected a {.cls tmle_spec} object, not {.cls {class(spec)}}.",
      call = call
    )
  }
  invisible(spec)
}

#' Configure Preprocessing Recipe for a TMLE Specification
#'
#' Attaches a base `recipes::recipe` to a [tmle_spec()]. During fitting,
#' `tidytmle` automatically clones and adapts this recipe per nuisance factor
#' (`Q`, `g`, `C`) so treatment \eqn{A} is never leaked into propensity models
#' and `id`/`censoring`/`weights` columns are excluded from predictors.
#'
#' @param spec A [tmle_spec()] object.
#' @param recipe A `recipes::recipe` object (or `NULL` to revert to automatic
#'   formula preprocessing).
#'
#' @return An updated `"tmle_spec"` object.
#'
#' @examples
#' df <- sim_causal_cohort(n = 120, seed = 1)
#' rec <- recipes::recipe(y ~ ., data = df) |>
#'   recipes::step_dummy(recipes::all_nominal_predictors())
#' spec <- tmle_spec(y, a, c(age, sex, lab_1)) |>
#'   set_preprocessing(rec)
#' spec
#'
#' @export
set_preprocessing <- function(spec, recipe) {
  check_is_tmle_spec(spec)
  if (!is.null(recipe) && !inherits(recipe, "recipe")) {
    cli::cli_abort("{.arg recipe} must be a {.cls recipe} object or {.code NULL}.")
  }
  spec$preprocessing <- recipe
  spec
}

validate_learners_arg <- function(learners, call = rlang::caller_env()) {
  if (is.null(learners) || identical(learners, "default")) {
    return("default")
  }
  valid <- inherits(learners, "model_spec") ||
    inherits(learners, "workflow") ||
    inherits(learners, "workflow_set") ||
    (is.list(learners) && length(learners) > 0L && all(vapply(learners, function(x) {
      inherits(x, "model_spec") || inherits(x, "workflow")
    }, logical(1L))))

  if (!valid) {
    cli::cli_abort(c(
      "{.arg learners} must be {.val default}, a {.cls parsnip} {.cls model_spec}, a {.cls workflow}, a {.cls workflow_set}, or a list of specs/workflows.",
      "x" = "Received object of class {.cls {class(learners)}}."
    ), call = call)
  }
  learners
}

#' Configure Nuisance Learners for Outcome, Treatment, or Censoring
#'
#' Polymorphically accepts a single `parsnip` model specification, a `workflow`,
#' a `workflowsets::workflow_set`, a list of specs/workflows, or `"default"`.
#' Passing multiple learners or a `workflow_set` automatically activates
#' cross-validated Super Learning within the exact same `rsample` folds.
#'
#' @param spec A [tmle_spec()] object.
#' @param learners `"default"`, a `parsnip` `model_spec`, a `workflows::workflow`,
#'   a `workflowsets::workflow_set`, or a list of `model_spec`/`workflow` objects.
#' @param metalearner Ensemble metalearning rule: `"auto"` (NNLS on simplex for
#'   regression, binomial log-loss minimization on simplex for classification),
#'   `"nnls"`, `"logloss"`, or `"best"` (discrete Super Learner selecting the
#'   single lowest CV-risk learner).
#' @param stratify Logical (for `set_outcome_learners()` with binary treatments);
#'   if `TRUE`, fits separate arm-stratified Super Learners \eqn{\bar{Q}(1, W)}
#'   on treated rows and \eqn{\bar{Q}(0, W)} on control rows (T-learner mode)
#'   instead of pooling treatment as a covariate. Defaults to `FALSE`.
#'
#' @return An updated `"tmle_spec"` object.
#'
#' @examples
#' spec <- tmle_spec(y, a, c(age, sex, lab_1)) |>
#'   set_outcome_learners(parsnip::linear_reg() |> parsnip::set_engine("glm")) |>
#'   set_treatment_learners(parsnip::logistic_reg() |> parsnip::set_engine("glm"))
#' spec
#'
#' @export
set_outcome_learners <- function(spec,
                                 learners,
                                 metalearner = c("auto", "nnls", "logloss", "best"),
                                 stratify = FALSE) {
  check_is_tmle_spec(spec)
  metalearner <- match.arg(metalearner)
  spec$learners$Q <- list(
    spec = validate_learners_arg(learners),
    metalearner = metalearner,
    stratify = isTRUE(stratify)
  )
  spec
}

#' @rdname set_outcome_learners
#' @export
set_treatment_learners <- function(spec,
                                   learners,
                                   metalearner = c("auto", "logloss", "nnls", "best")) {
  check_is_tmle_spec(spec)
  metalearner <- match.arg(metalearner)
  spec$learners$g <- list(
    spec = validate_learners_arg(learners),
    metalearner = metalearner
  )
  spec
}

#' @rdname set_outcome_learners
#' @export
set_censoring_learners <- function(spec,
                                   learners,
                                   metalearner = c("auto", "logloss", "nnls", "best")) {
  check_is_tmle_spec(spec)
  metalearner <- match.arg(metalearner)
  spec$learners$C <- list(
    spec = validate_learners_arg(learners),
    metalearner = metalearner
  )
  spec
}

#' Configure Cross-Fitting Resamples for a TMLE Specification
#'
#' @param spec A [tmle_spec()] object.
#' @param resamples An `rsample::rset` object (e.g., `vfold_cv()` or
#'   `group_vfold_cv()`), or `NULL` to auto-generate smart stratified/grouped
#'   resamples at fit time.
#' @param v_folds Integer number of cross-fitting folds when `resamples = NULL`.
#'   Defaults to `10L`.
#' @param strata Logical; whether to stratify auto-generated folds on treatment
#'   (and binary outcome when balanced). Defaults to `TRUE`.
#'
#' @return An updated `"tmle_spec"` object.
#'
#' @examples
#' spec <- tmle_spec(y, a, c(age, sex, lab_1)) |>
#'   set_resamples(v_folds = 5L)
#' spec
#'
#' @export
set_resamples <- function(spec, resamples = NULL, v_folds = 10L, strata = TRUE) {
  check_is_tmle_spec(spec)
  if (!is.null(resamples) && !inherits(resamples, "rset")) {
    cli::cli_abort("{.arg resamples} must be an {.cls rset} object from {.pkg rsample} or {.code NULL}.")
  }
  if (!is.numeric(v_folds) || length(v_folds) != 1L || v_folds < 2L) {
    cli::cli_abort("{.arg v_folds} must be an integer >= 2.")
  }
  spec$resamples <- list(
    obj = resamples,
    v_folds = as.integer(v_folds),
    strata = isTRUE(strata)
  )
  spec
}

#' Configure TMLE Targeting, Positivity Bounds, and Memory Options
#'
#' @details
#' Three built-in statistical guardrails govern the targeting step:
#' - **CV-TMLE (`cv_tmle = TRUE`)**: Uses out-of-fold nuisance predictions so
#'   standard errors remain valid even when using non-Donsker machine learning
#'   learners such as random forests or boosted trees (Zheng & van der Laan,
#'   2011).
#' - **Adaptive Positivity Truncation (`g_bounds = "adaptive"`)**: Prevents a
#'   handful of near-0 or near-1 propensity scores from inflating inverse weights
#'   and EIF variance, following the sample-size-aware \eqn{5 / (\sqrt{n}\ln n)}
#'   truncation rule (Petersen et al., 2012; Gruber et al., 2022).
#' - **Bounded Logistic Fluctuation (`fluctuation = "logistic"`)**: Scales
#'   continuous outcomes into `[0, 1]` and performs canonical logistic
#'   fluctuation so targeted counterfactual predictions stay strictly inside the
#'   observed support even under extreme confounding (Gruber & van der Laan,
#'   2010).
#'
#' @references
#' Zheng, W., & van der Laan, M. J. (2011). Cross-validated targeted minimum-loss-based
#' estimation. In *Targeted Learning* (pp. 459–474). Springer.
#'
#' Gruber, S., & van der Laan, M. J. (2010). A targeted maximum likelihood
#' estimator of a causal effect on a bounded continuous outcome. *The
#' International Journal of Biostatistics*, 6(1), Article 26.
#'
#' Gruber, S., Phillips, R. V., Lee, H., & van der Laan, M. J. (2022). Data-adaptive
#' selection of the propensity score truncation level for inverse-probability-weighted
#' and targeted maximum likelihood estimators of marginal point treatment effects.
#' *American Journal of Epidemiology*, 191(9), 1640–1651.
#'
#' Petersen, M. L., Porter, K. E., Gruber, S., Wang, Y., & van der Laan, M. J.
#' (2012). Diagnosing and responding to violations in the positivity assumption.
#' *Statistical Methods in Medical Research*, 21(1), 31–54.
#'
#' @param spec A [tmle_spec()] object.
#' @param cv_tmle Logical; whether to use Cross-Validated TMLE (out-of-fold
#'   nuisance predictions). Strongly recommended (`TRUE` by default) to preserve
#'   valid influence-function standard errors without Donsker assumptions.
#' @param g_bounds Positivity truncation rule for propensity scores: `"adaptive"`
#'   (sample-size-aware \eqn{5 / (\sqrt{n_{\text{eff}}} \ln n_{\text{eff}})} rule
#'   clamped to `[0.01, 0.15]`; Gruber et al., 2022), a single numeric threshold
#'   `b` giving `[b, 1 - b]`, or a length-2 numeric vector `c(low, high)` in
#'   `(0, 1)`.
#' @param fluctuation Fluctuation submodel: `"logistic"` (arm-specific
#'   weighted-intercept bounded canonical logistic fluctuation on scaled `[0, 1]`
#'   outcomes, equivalent to `target.gwt = TRUE` -- default),
#'   `"logistic_covariate"` (classical unweighted two-parameter Clever-Covariate
#'   logistic submodel matching `tmle::tmle(..., target.gwt = FALSE)`),
#'   `"linear"` (weighted linear fluctuation), or `"none"` (untargeted
#'   G-computation plug-in using initial ML predictions without a targeting
#'   fluctuation step; EIF score is not solved).
#' @param max_iter Maximum iterative TMLE steps (reserved for iterative updates).
#' @param tol Convergence tolerance reserved for iterative longitudinal updates
#'   (point-treatment canonical logistic updates converge in one closed-form GLM
#'   iteration using internal tolerance `1e-10`).
#' @param save_nuisance Memory policy for nuisance models after fitting:
#'   `"predictions"` (default; keeps only out-of-fold predictions, CV risk, and
#'   Super Learner weights for instant `retarget()` with minimal RAM), `"butcher"`
#'   (strips heavy internals via `butcher::butcher()`), or `"full"` (retains full
#'   fitted workflows).
#' @param parallel Logical; if `TRUE` and `mirai` daemons are active
#'   (`mirai::daemons_set()`), distributes `(fold x learner)` tasks via
#'   `mirai::mirai_map()`.
#'
#' @return An updated `"tmle_spec"` object.
#'
#' @examples
#' spec <- tmle_spec(y, a, c(age, sex, lab_1)) |>
#'   set_targeting(cv_tmle = TRUE, g_bounds = c(0.05, 0.95), fluctuation = "logistic")
#' spec
#'
#' @export
set_targeting <- function(spec,
                          cv_tmle = TRUE,
                          g_bounds = "adaptive",
                          fluctuation = c("logistic", "logistic_covariate", "linear", "none"),
                          max_iter = 100L,
                          tol = 1e-6,
                          save_nuisance = c("predictions", "butcher", "full"),
                          parallel = FALSE) {
  check_is_tmle_spec(spec)
  fluctuation <- match.arg(fluctuation)
  save_nuisance <- match.arg(save_nuisance)
  validate_g_bounds_arg(g_bounds)

  spec$targeting <- list(
    cv_tmle = isTRUE(cv_tmle),
    g_bounds = g_bounds,
    fluctuation = fluctuation,
    max_iter = as.integer(max_iter),
    tol = as.numeric(tol),
    save_nuisance = save_nuisance,
    parallel = isTRUE(parallel)
  )
  spec
}

validate_g_bounds_arg <- function(g_bounds, call = rlang::caller_env()) {
  if (identical(g_bounds, "adaptive")) {
    return(invisible(g_bounds))
  }
  if (is.numeric(g_bounds) && length(g_bounds) == 1L && !is.na(g_bounds) &&
      g_bounds > 0 && g_bounds < 0.5) {
    return(invisible(g_bounds))
  }
  if (is.numeric(g_bounds) && length(g_bounds) == 2L && !anyNA(g_bounds) &&
      g_bounds[1L] > 0 && g_bounds[1L] < g_bounds[2L] && g_bounds[2L] < 1) {
    return(invisible(g_bounds))
  }
  cli::cli_abort(c(
    "Invalid {.arg g_bounds} specification.",
    "i" = "Supply {.val adaptive}, a single scalar in (0, 0.5), or a pair {.code c(low, high)} with 0 < low < high < 1."
  ), call = call)
}
