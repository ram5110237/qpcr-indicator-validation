# ---------------------------------------------------------------------------
# marker_validation.R
# Diagnostic accuracy of molecular markers against a reference method
# (e.g. a qPCR/allele-specific assay vs Sanger or WGS genotype, or vs a
# bioassay phenotype), ROC analysis, agreement, and cross-validated
# performance of multi-marker models.
# ---------------------------------------------------------------------------

wilson_ci <- function(x, n, conf = 0.95) {
  if (n == 0) return(c(NA, NA))
  z <- stats::qnorm(1 - (1 - conf) / 2); p <- x / n
  centre <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  c(centre - half, centre + half)
}

#' Sensitivity, specificity, predictive values, accuracy and Cohen's kappa
diagnostic_metrics <- function(test, reference) {
  test <- as.logical(test); reference <- as.logical(reference)
  ok <- !is.na(test) & !is.na(reference)
  t <- test[ok]; r <- reference[ok]
  tp <- sum(t & r); tn <- sum(!t & !r); fp <- sum(t & !r); fn <- sum(!t & r); n <- length(t)
  row <- function(name, x, m) { ci <- wilson_ci(x, m); data.frame(metric = name, estimate = if (m) x / m else NA,
                                                                    ci_low = ci[1], ci_high = ci[2], n = m) }
  po <- (tp + tn) / n
  pe <- ((tp + fp) * (tp + fn) + (fn + tn) * (fp + tn)) / n^2
  kappa <- if (pe < 1) (po - pe) / (1 - pe) else NA
  se_k <- sqrt(po * (1 - po) / (n * (1 - pe)^2))
  rbind(
    row("sensitivity", tp, tp + fn), row("specificity", tn, tn + fp),
    row("PPV", tp, tp + fp), row("NPV", tn, tn + fn), row("accuracy", tp + tn, n),
    data.frame(metric = "cohen_kappa", estimate = kappa, ci_low = kappa - 1.96 * se_k,
               ci_high = kappa + 1.96 * se_k, n = n)
  )
}

#' ROC curve, AUC with DeLong CI, and the Youden-optimal threshold
roc_marker <- function(score, truth, direction = "auto") {
  r <- pROC::roc(response = as.integer(as.logical(truth)), predictor = score,
                 levels = c(0, 1), direction = if (direction == "auto") "auto" else direction, quiet = TRUE)
  ci <- as.numeric(pROC::ci.auc(r, method = "delong"))
  best <- pROC::coords(r, "best", best.method = "youden",
                       ret = c("threshold", "sensitivity", "specificity"), transpose = FALSE)[1, ]
  list(roc = r,
       summary = data.frame(auc = ci[2], auc_ci_low = ci[1], auc_ci_high = ci[3],
                            youden_threshold = best$threshold, sens_at_threshold = best$sensitivity,
                            spec_at_threshold = best$specificity, direction = r$direction))
}

#' Repeated stratified k-fold cross-validated AUC of a logistic model
cv_auc <- function(data, predictors, outcome, k = 5, repeats = 20, seed = 1) {
  set.seed(seed)
  y <- as.integer(as.logical(data[[outcome]]))
  f <- stats::as.formula(paste(outcome, "~", paste(predictors, collapse = " + ")))
  d <- data; d[[outcome]] <- y
  aucs <- replicate(repeats, {
    folds <- integer(nrow(d))
    for (cl in 0:1) { idx <- which(y == cl); folds[idx] <- sample(rep_len(1:k, length(idx))) }
    pred <- numeric(nrow(d))
    for (i in 1:k) {
      fit <- suppressWarnings(stats::glm(f, family = stats::binomial, data = d[folds != i, ]))
      pred[folds == i] <- stats::predict(fit, d[folds == i, ], type = "response")
    }
    as.numeric(pROC::auc(pROC::roc(y, pred, levels = c(0, 1), direction = "<", quiet = TRUE)))
  })
  data.frame(model = paste(predictors, collapse = " + "), cv_auc_mean = mean(aucs),
             cv_auc_sd = stats::sd(aucs), k = k, repeats = repeats)
}

#' Rank individual markers and compare with the combined model
marker_panel_evaluation <- function(data, markers, outcome, k = 5, repeats = 20) {
  single <- do.call(rbind, lapply(markers, function(m) {
    r <- roc_marker(data[[m]], data[[outcome]])$summary
    cbind(marker = m, r, cv_auc(data, m, outcome, k, repeats)[, c("cv_auc_mean", "cv_auc_sd")])
  }))
  single <- single[order(-single$auc), ]
  combined <- cv_auc(data, markers, outcome, k, repeats)
  list(single = single, combined = combined)
}
