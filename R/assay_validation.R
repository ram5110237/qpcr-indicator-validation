# ---------------------------------------------------------------------------
# assay_validation.R
# qPCR assay performance: standard curve and efficiency, limit of detection,
# limit of quantification, repeatability / intermediate precision, and
# efficiency-corrected relative quantification with multiple reference genes.
# Terminology follows the MIQE guidelines (Bustin et al. 2009).
# ---------------------------------------------------------------------------

#' Standard curve: Cq ~ log10(quantity) over detected reactions
standard_curve <- function(quantity, cq) {
  ok <- !is.na(cq) & quantity > 0
  if (length(unique(quantity[ok])) < 3) stop("Need at least 3 dilution levels with detected Cq")
  fit <- stats::lm(cq[ok] ~ log10(quantity[ok]))
  slope <- unname(stats::coef(fit)[2])
  ci <- stats::confint(fit)[2, ]
  eff <- 10^(-1 / slope) - 1
  data.frame(
    slope = slope, intercept = unname(stats::coef(fit)[1]),
    r_squared = summary(fit)$r.squared,
    efficiency_pct = 100 * eff,
    efficiency_ci_low = 100 * (10^(-1 / ci[1]) - 1),
    efficiency_ci_high = 100 * (10^(-1 / ci[2]) - 1),
    amplification_factor = 1 + eff,
    n_points = sum(ok),
    passes_miqe = (eff >= 0.9 & eff <= 1.1) & summary(fit)$r.squared >= 0.98
  )
}

#' Limit of detection: concentration detected with 95% probability (probit)
#'
#' Fitted on log10 concentration; the CI comes from the delta method on the
#' log scale (MASS::dose.p) and is back-transformed.
lod_probit <- function(quantity, detected, p = 0.95) {
  d <- data.frame(lq = log10(quantity), det = as.integer(detected))
  if (all(d$det == 1)) {
    warning("All reactions detected: LOD is below the lowest tested level")
    return(data.frame(lod = min(quantity), lod_ci_low = NA, lod_ci_high = NA, method = "all_detected"))
  }
  fit <- suppressWarnings(stats::glm(det ~ lq, family = stats::binomial("probit"), data = d))
  dp <- MASS::dose.p(fit, p = p)
  est <- as.numeric(dp); se <- as.numeric(attr(dp, "SE"))
  data.frame(lod = 10^est, lod_ci_low = 10^(est - 1.96 * se), lod_ci_high = 10^(est + 1.96 * se),
             method = paste0("probit_", 100 * p))
}

#' Detection rate per concentration level (useful alongside the probit LOD)
detection_table <- function(quantity, detected) {
  a <- stats::aggregate(detected ~ quantity, FUN = function(x) c(n = length(x), detected = sum(x)))
  data.frame(quantity = a$quantity, n = a$detected[, "n"], detected = a$detected[, "detected"],
             rate = a$detected[, "detected"] / a$detected[, "n"])
}

#' Limit of quantification: lowest level whose back-calculated quantities have
#' CV <= max_cv (MIQE-style default 35%), with all reactions detected.
loq_cv <- function(quantity, cq, curve, max_cv = 35) {
  back <- 10^((cq - curve$intercept) / curve$slope)
  lev <- sort(unique(quantity))
  tab <- do.call(rbind, lapply(lev, function(q) {
    b <- back[quantity == q]
    data.frame(quantity = q, n = length(b), all_detected = !anyNA(b),
               cv_pct = if (sum(!is.na(b)) >= 2) 100 * stats::sd(b, na.rm = TRUE) / mean(b, na.rm = TRUE) else NA)
  }))
  tab$quantifiable <- tab$all_detected & !is.na(tab$cv_pct) & tab$cv_pct <= max_cv
  # LOQ = lowest level from which all higher levels are also quantifiable
  ok <- rev(cumprod(rev(tab$quantifiable)) == 1)
  list(table = tab, loq = if (any(ok)) min(tab$quantity[ok]) else NA)
}

#' Precision from a dilution series run on several plates/days
#'
#' Mixed model per quantity level: Cq ~ (1 | run). Repeatability = within-run SD,
#' intermediate precision = sqrt(within + between run variance).
precision_components <- function(df, min_detect_frac = 1) {
  out <- lapply(split(df, df$quantity), function(d) {
    if (mean(!is.na(d$cq)) < min_detect_frac || length(unique(d$run)) < 2) return(NULL)
    d <- d[!is.na(d$cq), ]
    fit <- suppressMessages(lme4::lmer(cq ~ 1 + (1 | run), data = d, REML = TRUE))
    vc <- as.data.frame(lme4::VarCorr(fit))
    s_run <- vc$sdcor[vc$grp == "run"]; s_res <- vc$sdcor[vc$grp == "Residual"]
    data.frame(quantity = d$quantity[1], mean_cq = mean(d$cq), sd_repeatability = s_res,
               sd_between_run = s_run, sd_intermediate = sqrt(s_res^2 + s_run^2),
               cv_repeatability_pct = 100 * s_res / mean(d$cq))
  })
  do.call(rbind, out)
}

#' Efficiency-corrected relative expression (Pfaffl 2001) with several
#' reference genes combined by geometric mean (Vandesompele et al. 2002).
#'
#' @param target_cq,ref_cq Cq for target and a matrix of reference genes
#' @param group factor with a control level; ratios are relative to its mean
#' @param e_target,e_ref amplification factors (2 = 100%)
relative_expression <- function(target_cq, ref_cq, group, control, e_target = 2, e_ref = NULL) {
  ref_cq <- as.matrix(ref_cq)
  if (is.null(e_ref)) e_ref <- rep(2, ncol(ref_cq))
  group <- stats::relevel(factor(group), ref = control)
  ctl <- group == control
  d_target <- (mean(target_cq[ctl]) - target_cq) * log2(e_target)
  d_ref <- sweep(sweep(-ref_cq, 2, colMeans(ref_cq[ctl, , drop = FALSE]), "+"), 2, log2(e_ref), "*")
  log2_ratio <- d_target - rowMeans(d_ref)
  fit <- stats::lm(log2_ratio ~ group)
  est <- summary(fit)$coefficients[-1, , drop = FALSE]
  ci <- stats::confint(fit)[-1, , drop = FALSE]
  list(per_sample = data.frame(group = group, log2_ratio = log2_ratio),
       contrasts = data.frame(comparison = paste(sub("^group", "", rownames(est)), "vs", control),
                              log2_fold_change = est[, 1], fold_change = 2^est[, 1],
                              ci_low = 2^ci[, 1], ci_high = 2^ci[, 2], pvalue = est[, 4],
                              row.names = NULL))
}
