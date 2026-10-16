# ---------------------------------------------------------------------------
# dose_response.R
# Fungicide sensitivity: four-parameter log-logistic fits, EC50 with CI,
# resistance factors relative to a reference (baseline) isolate.
#
# Model (same parameterisation as drc::LL.4):
#   y = c + (d - c) / (1 + exp(b * (log(x) - log(e))))
# with e = EC50. log(e) is estimated directly so the CI is symmetric on the
# log scale and always positive when back-transformed.
# ---------------------------------------------------------------------------

ll4 <- function(x, b, c, d, loge) c + (d - c) / (1 + exp(b * (log(x) - loge)))

#' Fit LL.4 to one isolate. Zero doses (untreated controls) are allowed.
fit_ll4 <- function(dose, response, fix_lower = NULL) {
  pos <- dose > 0
  start <- list(b = 1, d = max(response), loge = log(stats::median(dose[pos])))
  if (is.null(fix_lower)) start$c <- min(response)
  frm <- if (is.null(fix_lower)) response ~ ll4(dose, b, c, d, loge)
         else response ~ ll4(dose, b, fix_lower, d, loge)
  # try a few slope/EC50 starts; nls is sensitive to starting values
  grid <- expand.grid(b = c(1, 2, 0.5), le = start$loge + c(0, -1, 1))
  fit <- NULL
  for (i in seq_len(nrow(grid))) {
    s <- start; s$b <- grid$b[i]; s$loge <- grid$le[i]
    fit <- tryCatch(stats::nls(frm, start = s, control = stats::nls.control(maxiter = 200, warnOnly = FALSE)),
                    error = function(e) NULL)
    if (!is.null(fit)) break
  }
  if (is.null(fit)) stop("LL.4 fit did not converge")
  fit
}

ec50_summary <- function(fit, level = 0.95) {
  co <- summary(fit)$coefficients
  le <- co["loge", "Estimate"]; se <- co["loge", "Std. Error"]
  z <- stats::qt(1 - (1 - level) / 2, df = summary(fit)$df[2])
  data.frame(ec50 = exp(le), ec50_ci_low = exp(le - z * se), ec50_ci_high = exp(le + z * se),
             log_ec50 = le, se_log_ec50 = se, slope_b = co["b", "Estimate"],
             upper_d = co["d", "Estimate"],
             lower_c = if ("c" %in% rownames(co)) co["c", "Estimate"] else NA,
             residual_se = summary(fit)$sigma)
}

#' Fit every isolate, then compute resistance factors vs `reference`
#'
#' RF = EC50_isolate / EC50_reference; its CI uses var(log RF) = se1^2 + se2^2
#' (independent fits). Isolates are labelled by RF class (thresholds editable).
fit_isolates <- function(df, reference, rf_breaks = c(0, 5, 10, 50, Inf),
                         rf_labels = c("sensitive", "low resistance", "moderate resistance", "high resistance"),
                         fix_lower = NULL) {
  res <- do.call(rbind, lapply(split(df, df$isolate), function(d) {
    s <- tryCatch(ec50_summary(fit_ll4(d$dose, d$response, fix_lower)),
                  error = function(e) data.frame(ec50 = NA, ec50_ci_low = NA, ec50_ci_high = NA,
                                                 log_ec50 = NA, se_log_ec50 = NA, slope_b = NA,
                                                 upper_d = NA, lower_c = NA, residual_se = NA))
    cbind(isolate = d$isolate[1], s)
  }))
  if (!reference %in% res$isolate) stop("Reference isolate not found: ", reference)
  ref <- res[res$isolate == reference, ]
  res$rf <- res$ec50 / ref$ec50
  se_rf <- sqrt(res$se_log_ec50^2 + ref$se_log_ec50^2)
  res$rf_ci_low <- exp(log(res$rf) - 1.96 * se_rf)
  res$rf_ci_high <- exp(log(res$rf) + 1.96 * se_rf)
  res$rf_ci_low[res$isolate == reference] <- res$rf_ci_high[res$isolate == reference] <- 1
  res$sensitivity_class <- as.character(cut(res$rf, rf_breaks, labels = rf_labels, right = FALSE))
  rownames(res) <- NULL
  res[order(res$rf), ]
}

#' Predicted curves for plotting
predict_curves <- function(df, fits_table, n = 100) {
  do.call(rbind, lapply(fits_table$isolate[!is.na(fits_table$ec50)], function(i) {
    f <- fits_table[fits_table$isolate == i, ]
    pos <- df$dose[df$isolate == i & df$dose > 0]
    x <- exp(seq(log(min(pos)), log(max(pos)), length.out = n))
    c0 <- if (is.na(f$lower_c)) 0 else f$lower_c
    data.frame(isolate = i, dose = x, response = ll4(x, f$slope_b, c0, f$upper_d, f$log_ec50))
  }))
}
