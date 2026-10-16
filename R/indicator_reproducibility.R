# ---------------------------------------------------------------------------
# indicator_reproducibility.R
# How reproducible is a quantitative indicator (marker expression, disease
# score, pathogen load) across environments, sites or years?
#
# Linear mixed model with crossed random effects:
#   y ~ 1 + (1 | entry) + (1 | env) + (1 | entry:env)
# "entry" can be a genotype, cultivar, isolate or sample source.
# ---------------------------------------------------------------------------

variance_components <- function(fit) {
  vc <- as.data.frame(lme4::VarCorr(fit))
  stats::setNames(vc$vcov, vc$grp)
}

#' Fit the model and derive repeatability / broad-sense heritability on an
#' entry-mean basis: H2 = Vg / (Vg + Vge/e + Ve/(e*r))
fit_indicator_model <- function(df, response = "value", entry = "entry", env = "env") {
  d <- data.frame(y = df[[response]], entry = factor(df[[entry]]), env = factor(df[[env]]))
  if (nlevels(d$env) < 2) stop("Need at least two environments")
  frm <- if (any(table(d$entry, d$env) > 1)) y ~ 1 + (1 | entry) + (1 | env) + (1 | entry:env)
         else y ~ 1 + (1 | entry) + (1 | env)
  fit <- suppressMessages(lme4::lmer(frm, data = d, REML = TRUE))
  list(fit = fit, data = d, summary = indicator_stats(fit, d))
}

indicator_stats <- function(fit, d) {
  v <- variance_components(fit)
  vg <- v[["entry"]]; vge <- if ("entry:env" %in% names(v)) v[["entry:env"]] else 0
  ve <- v[["Residual"]]; venv <- v[["env"]]
  n_env <- nlevels(d$env)
  n_rep <- mean(table(d$entry, d$env))
  h2 <- vg / (vg + vge / n_env + ve / (n_env * n_rep))
  icc <- vg / (vg + vge + ve)                      # single-observation reliability
  data.frame(var_entry = vg, var_env = venv, var_entry_env = vge, var_residual = ve,
             n_env = n_env, mean_reps = n_rep, h2_entry_mean = h2, icc_single = icc,
             ge_ratio = if (vg > 0) vge / vg else NA)
}

#' Parametric bootstrap CI for H2 and ICC
bootstrap_indicator <- function(model, nsim = 200, seed = 1) {
  fn <- function(f) {
    s <- indicator_stats(f, model$data)
    c(h2 = s$h2_entry_mean, icc = s$icc_single)
  }
  b <- suppressMessages(suppressWarnings(lme4::bootMer(model$fit, fn, nsim = nsim, seed = seed)))
  data.frame(statistic = c("h2_entry_mean", "icc_single"),
             estimate = b$t0,
             ci_low = apply(b$t, 2, stats::quantile, 0.025, na.rm = TRUE),
             ci_high = apply(b$t, 2, stats::quantile, 0.975, na.rm = TRUE), row.names = NULL)
}

#' Entry BLUPs (shrunken entry effects) with conditional SE
entry_blups <- function(model) {
  re <- lme4::ranef(model$fit, condVar = TRUE)$entry
  pv <- attr(re, "postVar")
  data.frame(entry = rownames(re), blup = re[, 1], se = sqrt(pv[1, 1, ]),
             predicted_mean = lme4::fixef(model$fit)[1] + re[, 1], row.names = NULL)[order(-re[, 1]), ]
}

#' Finlay-Wilkinson stability: regress each entry's env mean on the env index
finlay_wilkinson <- function(df, response = "value", entry = "entry", env = "env") {
  m <- stats::aggregate(df[[response]], list(entry = df[[entry]], env = df[[env]]), mean)
  names(m)[3] <- "y"
  idx <- stats::aggregate(y ~ env, m, mean); idx$env_index <- idx$y - mean(idx$y)
  m <- merge(m, idx[, c("env", "env_index")])
  do.call(rbind, lapply(split(m, m$entry), function(d) {
    f <- stats::lm(y ~ env_index, d)
    data.frame(entry = d$entry[1], mean = mean(d$y), fw_slope = unname(stats::coef(f)[2]),
               deviation_mse = summary(f)$sigma^2)
  }))
}
