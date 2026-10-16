repo <- normalizePath(file.path(testthat::test_path(), "..", ".."))
for (f in list.files(file.path(repo, "R"), pattern = "\\.R$", full.names = TRUE)) source(f)

sim_ct <- function(n = 60, sds = c(0.1, 0.2, 0.4, 0.8), shift = NULL, seed = 1) {
  set.seed(seed)
  inp <- rnorm(n, 0, 1)
  ct <- sapply(sds, function(s) 22 + inp + rnorm(n, 0, s))
  colnames(ct) <- paste0("g", seq_along(sds))
  ct
}

# ---------------------------------------------------------------- reference genes
test_that("geNorm M equals the mean SD of pairwise log-ratios", {
  ct <- sim_ct()
  gn <- genorm(ct)
  m1 <- mean(sapply(2:4, function(k) sd(ct[, 1] - ct[, k])))
  first_step <- gn$M_steps[gn$M_steps$n_genes == 4, ]
  expect_equal(first_step$M[first_step$gene == "g1"], m1)
  expect_equal(gn$ranking$gene[4], "g4")              # noisiest excluded first
  expect_true(all(c("g1", "g2") %in% gn$ranking$gene[1:2]))
})

test_that("efficiency below 100% shrinks log2 ratios proportionally", {
  ct <- sim_ct()
  expect_equal(log2_rq(ct, rep(1.9, 4)), -ct * log2(1.9))
})

test_that("NormFinder recovers true intragroup variances and penalises a group shift", {
  set.seed(2)
  n <- 400; g <- rep(c("a", "b"), each = n / 2)
  ct <- sim_ct(n, sds = c(0.1, 0.2, 0.3, 0.5), seed = 3)
  nf1 <- normfinder(ct)
  expect_equal(nf1$gene, paste0("g", 1:4))
  expect_equal(sqrt(nf1$intragroup_var[match(paste0("g", 1:4), nf1$gene)]),
               c(0.1, 0.2, 0.3, 0.5), tolerance = 0.25)
  ct[g == "b", "g1"] <- ct[g == "b", "g1"] + 1       # most precise gene, but condition-dependent
  nf2 <- normfinder(ct, g)
  expect_equal(nf2$gene[4], "g1")
})

test_that("BestKeeper and delta-Ct rank the noisiest gene last", {
  ct <- sim_ct()
  expect_equal(tail(bestkeeper(ct)$gene, 1), "g4")
  expect_equal(tail(delta_ct(ct)$gene, 1), "g4")
  expect_equal(consensus_ranking(ct)$consensus$gene[1], "g1")
})

test_that("missing Cq values are rejected", {
  ct <- sim_ct(); ct[1, 1] <- NA
  expect_error(genorm(ct), "missing")
})

# ---------------------------------------------------------------- assay validation
test_that("standard curve recovers slope and efficiency", {
  q <- rep(10^(6:2), each = 3)
  cq <- 38 - 3.32 * log10(q) + rnorm(length(q), 0, 0.01)
  sc <- standard_curve(q, cq)
  expect_equal(sc$slope, -3.32, tolerance = 0.01)
  expect_equal(sc$efficiency_pct, 100 * (10^(1 / 3.32) - 1), tolerance = 0.5)
  expect_true(sc$passes_miqe)
})

test_that("probit LOD is close to the Poisson theoretical value", {
  set.seed(4)
  q <- rep(c(0.5, 1, 2, 3, 5, 8, 10, 20), each = 200)
  det <- rpois(length(q), q) > 0
  lod <- lod_probit(q, det)
  # P(detect) = 1 - exp(-q) = 0.95  ->  q = 3.0 copies
  expect_gt(lod$lod, 2); expect_lt(lod$lod, 5.5)
  expect_lt(lod$lod_ci_low, lod$lod); expect_gt(lod$lod_ci_high, lod$lod)
})

test_that("LOQ requires all higher levels to be quantifiable", {
  curve <- data.frame(intercept = 38, slope = -3.32)
  q <- rep(c(1000, 100, 10, 1), each = 4)
  cq <- 38 - 3.32 * log10(q) + c(rep(0.05, 8) * c(1, -1), c(0.5, -0.6, 0.4, -0.3), c(1.5, NA, -1, 0.5))
  l <- loq_cv(q, cq, curve, 35)
  expect_equal(l$loq, 100)
})

test_that("precision components separate run and residual variance", {
  set.seed(5)
  d <- expand.grid(rep = 1:20, run = paste0("r", 1:8), quantity = 1000)
  d$cq <- 28 + rnorm(8, 0, 0.3)[as.integer(factor(d$run))] + rnorm(nrow(d), 0, 0.1)
  p <- precision_components(d)
  expect_equal(p$sd_repeatability, 0.1, tolerance = 0.2)
  expect_gt(p$sd_between_run, p$sd_repeatability)
})

test_that("relative expression recovers a planted fold change with two references", {
  set.seed(6)
  grp <- rep(c("ctl", "trt"), each = 10); inp <- rnorm(20, 0, 1)
  refs <- cbind(r1 = 20 + inp + rnorm(20, 0, 0.1), r2 = 24 + inp + rnorm(20, 0, 0.1))
  target <- 27 + inp - 2 * (grp == "trt") + rnorm(20, 0, 0.1)
  re <- relative_expression(target, refs, grp, "ctl")
  expect_equal(re$contrasts$fold_change, 4, tolerance = 0.1)
})

# ---------------------------------------------------------------- marker validation
test_that("diagnostic metrics match a hand-computed 2x2 table", {
  test <- c(rep(TRUE, 45), rep(FALSE, 5), rep(TRUE, 10), rep(FALSE, 40))
  ref  <- c(rep(TRUE, 50), rep(FALSE, 50))
  dm <- diagnostic_metrics(test, ref)
  expect_equal(dm$estimate[dm$metric == "sensitivity"], 0.9)
  expect_equal(dm$estimate[dm$metric == "specificity"], 0.8)
  expect_equal(dm$estimate[dm$metric == "cohen_kappa"], 0.7)
  ci <- dm[dm$metric == "sensitivity", c("ci_low", "ci_high")]
  expect_true(ci$ci_low < 0.9 && ci$ci_high > 0.9 && ci$ci_high <= 1)
})

test_that("ROC handles markers where low scores indicate the positive class", {
  set.seed(7)
  y <- rep(c(TRUE, FALSE), each = 50)
  s <- ifelse(y, rnorm(100, -2), rnorm(100, 2))
  r <- roc_marker(s, y)
  expect_gt(r$summary$auc, 0.9)
  expect_equal(r$summary$direction, ">")
})

test_that("cross-validated AUC is near 0.5 for a pure-noise marker", {
  set.seed(8)
  d <- data.frame(y = rep(c(TRUE, FALSE), 50), noise = rnorm(100))
  cv <- cv_auc(d, "noise", "y", k = 5, repeats = 10)
  expect_lt(abs(cv$cv_auc_mean - 0.5), 0.15)
})

# ---------------------------------------------------------------- dose response
test_that("LL.4 fit recovers EC50 and resistance factors", {
  set.seed(9)
  doses <- c(0, 0.01, 0.03, 0.1, 0.3, 1, 3, 10)
  mk <- function(iso, e) data.frame(isolate = iso, dose = rep(doses, 3),
                                    response = 100 / (1 + exp(1.5 * (log(pmax(rep(doses, 3), 1e-12)) - log(e)))) +
                                      rnorm(24, 0, 2))
  df <- rbind(mk("ref", 0.05), mk("res", 1))
  fits <- fit_isolates(df, "ref", fix_lower = 0)
  expect_equal(fits$ec50[fits$isolate == "ref"], 0.05, tolerance = 0.15)
  expect_equal(fits$rf[fits$isolate == "res"], 20, tolerance = 0.2)
  expect_equal(fits$sensitivity_class[fits$isolate == "res"], "moderate resistance")
  expect_true(all(fits$ec50_ci_low < fits$ec50 & fits$ec50 < fits$ec50_ci_high))
})

# ---------------------------------------------------------------- multi-environment
test_that("variance components and H2 are recovered from a balanced trial", {
  set.seed(10)
  ent <- sprintf("g%02d", 1:60); env <- paste0("e", 1:5)
  d <- expand.grid(rep = 1:3, env = env, entry = ent, stringsAsFactors = FALSE)
  ge <- rnorm(300, 0, sqrt(0.25)); names(ge) <- paste(rep(ent, each = 5), env)
  d$value <- 10 + rnorm(60, 0, 1)[match(d$entry, ent)] + rnorm(5, 0, 1)[match(d$env, env)] +
    ge[paste(d$entry, d$env)] + rnorm(nrow(d), 0, sqrt(0.5))
  m <- fit_indicator_model(d)
  s <- m$summary
  expect_equal(s$var_entry, 1, tolerance = 0.4)
  expect_equal(s$var_residual, 0.5, tolerance = 0.2)
  h2_true <- 1 / (1 + 0.25 / 5 + 0.5 / 15)
  expect_equal(s$h2_entry_mean, h2_true, tolerance = 0.05)
  fw <- finlay_wilkinson(d)
  expect_equal(mean(fw$fw_slope), 1, tolerance = 1e-8)   # FW slopes average to 1 by construction
})
