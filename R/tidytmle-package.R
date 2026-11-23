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
  ".estimate", ".estimator", ".metric", ".pred", "arm", "cate_val",
  "conf.high", "conf.low", "covariate", "cum_gain", "curve", "density",
  "eif", "estimand", "estimate", "fold", "frac", "g", "gain", "group",
  "h_mag", "is_outlier", "learner", "leverage", "metric_type", "nuisance",
  "obs_id", "pct_truncated", "random_gain", "smd", "term", "threshold",
  "val", "value", "weight", "x_init", "y_star"
))
