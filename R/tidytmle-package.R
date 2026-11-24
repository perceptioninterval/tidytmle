#' @keywords internal
"_PACKAGE"

#' @importFrom rlang enquo quo_is_null quo_is_missing quo_text eval_tidy as_label sym caller_env :=
#' @importFrom generics fit tidy glance augment
#' @importFrom ggplot2 autoplot
#' @importFrom tune collect_metrics
NULL

#' @export
generics::fit

#' @export
generics::tidy

#' @export
generics::glance

#' @export
generics::augment

#' @export
ggplot2::autoplot

#' @export
tune::collect_metrics

utils::globalVariables(c(
  ".estimate", ".estimator", ".metric", ".pred", "arm", "benefit_group",
  "cate_val", "conf.high", "conf.low", "covariate", "cum_gain", "curve",
  "density", "direct_lbl", "eif", "estimand", "estimate", "flag", "fold",
  "frac", "g", "gain", "group", "h_mag", "is_obs_contrast", "is_outlier", "label_human",
  "learner", "leverage", "mean", "metric_type", "null_ref", "nuisance",
  "nuisance_order", "obs_id", "pct_truncated", "plugin_est", "random_gain",
  "role_fill", "row_panel", "scale_group", "smd", "smd_raw", "smd_wt",
  "term", "threshold", "val", "val_lbl", "value", "weight", "x_init", "y_star"
))
