## =============================================================================
## slm_functions.R
## Exact covariate-adjusted linear regression from stratified summary statistics
##
## Input: one row per stratum, with the stratum-defining covariates, the stratum
## size n_i, the stratum sample mean ybar_i, and the stratum sample SD s_i.
## =============================================================================

## -----------------------------------------------------------------------------
## Internal: validate and prepare summary data
## -----------------------------------------------------------------------------
.slm_prepare <- function(formula, data, n, mean, sd) {
  for (v in c(n, mean, sd)) {
    if (!v %in% names(data)) stop("Column '", v, "' not found in data.")
  }
  data = data[!is.na(data[[n]]) & data[[n]] > 0, , drop = FALSE]
  if (any(is.na(data[[mean]]))) stop("Missing stratum means for non-empty strata.")

  ni = data[[n]]
  s  = data[[sd]]
  ## A stratum with n_i = 1 has an undefined SD but contributes (n_i - 1) S_i^2 = 0.
  s[ni == 1] = 0
  if (any(is.na(s))) stop("Missing SD for a stratum with n_i >= 2.")

  ## Rows with identical covariate patterns are allowed; they are treated as
  ## separate strata (the estimator is unaffected; M counts rows).
  X = model.matrix(formula, data)
  if (qr(X)$rank < ncol(X)) {
    stop("Design matrix is rank deficient: some coefficients are not identified ",
         "from the observed strata (check for empty covariate levels).")
  }
  list(data = data, X = X, ni = ni, ybar = data[[mean]], s = s)
}

## -----------------------------------------------------------------------------
## slm(): summary-data linear model
##
## formula : one-sided formula for the covariates, e.g. ~ covid + sex + race + age
## data    : data frame, one row per stratum
## n, mean, sd : column names for stratum size, mean, and SD
##
## Returns estimates identical to lm() fitted to the (unobserved) individual
## data, together with model-based standard errors, which are also exact individual-level quantities.
## -----------------------------------------------------------------------------
slm <- function(formula, data, n = "n", mean = "mean", sd = "sd") {
  p    = .slm_prepare(formula, data, n, mean, sd)
  X    = p$X; ni <- p$ni; ybar <- p$ybar; s <- p$s
  M    = nrow(X); K <- ncol(X); N <- sum(ni)
  if (N - K < 1) stop("N - K must be positive.")

  XtWX  = crossprod(X, X * ni)
  V     = solve(XtWX)                                  ## (X'WX)^{-1} (model-based cov)
  beta  = drop(V %*% crossprod(X, ni * ybar))
  r     = drop(ybar - X %*% beta)                      ## stratum-mean residuals

  SSE0  = sum((ni - 1) * s^2)                          ## within-stratum SS
  SSE1  = sum(ni * r^2)                                ## weighted between-stratum SS
  df    = N - K
  sigma2 = (SSE0 + SSE1) / df

  # Model-based covariance matrix
  vc_model = sigma2 * V

  dimnames(vc_model) = list(colnames(X), colnames(X))

  structure(list(
    coefficients = beta, vcov = vc_model, 
    sigma2 = sigma2, df = df, N = N, M = M, K = K,
    SSE0 = SSE0, SSE1 = SSE1, residuals_mean = r,
    X = X, n = ni, ybar = ybar, s = s, formula = formula
  ), class = "slm")
}

## Coefficient table
coef_table <- function(fit, level = 0.95) {
  vc = fit$vcov
  est = fit$coefficients
  se  = sqrt(diag(vc))
  tval = est / se
  pval = 2 * pt(-abs(tval), df = fit$df)
  q   = qt(1 - (1 - level) / 2, df = fit$df)
  data.frame(term = names(est), estimate = est, se = se, t = tval, p = pval,
             lower = est - q * se, upper = est + q * se,
             row.names = NULL, check.names = FALSE)
}

print.slm <- function(x, digits = 4, ...) {
  cat("Summary-data linear model (sLM)\n")
  cat(sprintf("  Strata M = %d, patients N = %d, parameters K = %d, df = %d\n",  x$M, x$N, x$K, x$df))
  cat(sprintf("  sigma-hat = %.4f\n\n", sqrt(x$sigma2)))
  print(coef_table(x), digits = digits)
  invisible(x)
}

## -----------------------------------------------------------------------------
## General linear hypothesis H0: A' beta = b0  (A is K x L)
## -----------------------------------------------------------------------------
## return: F-stat, df1/df2; p-value
slm_ftest <- function(fit, A, b0 = rep(0, ncol(A))) {
  A = as.matrix(A)
  if (nrow(A) != fit$K) stop("A must have K = ", fit$K, " rows.")
  L  = ncol(A)
  d  = drop(crossprod(A, fit$coefficients) - b0)
  vc = fit$vcov
  Fs = drop(t(d) %*% solve(t(A) %*% vc %*% A) %*% d) / L
  data.frame(F = Fs, df1 = L, df2 = fit$df, p = pf(Fs, L, fit$df, lower.tail = FALSE))
}

## Convenience: test that a set of named coefficients are all zero
slm_test_terms <- function(fit, terms) {
  idx = match(terms, names(fit$coefficients))
  if (anyNA(idx)) stop("Unknown term(s): ", paste(terms[is.na(idx)], collapse = ", "))
  A = diag(fit$K)[, idx, drop = FALSE]
  slm_ftest(fit, A)
}

## -----------------------------------------------------------------------------
## Diagnostics available from summary data
## -----------------------------------------------------------------------------

## Lack-of-fit test: compares between-stratum residual variation (SSE1) with pure error (SSE0). Requires M > K and N > M.
slm_lack_of_fit <- function(fit) {
  if (fit$M <= fit$K) return(data.frame(F = NA, df1 = 0, df2 = fit$N - fit$M, p = NA))
  df1 = fit$M - fit$K; df2 = fit$N - fit$M
  Fs  = (fit$SSE1 / df1) / (fit$SSE0 / df2)
  data.frame(F = Fs, df1 = df1, df2 = df2, p = pf(Fs, df1, df2, lower.tail = FALSE))
}

## Bartlett's test of equal within-stratum variances, from (n_i, s_i).
## Note: sensitive to non-normality; with skewed outcomes, a
## rejection may reflect skewness as much as heteroskedasticity.
slm_bartlett <- function(fit) {
  keep = (fit$n >= 2) & (fit$s > 0)
  v  = fit$n[keep] - 1; s2 = fit$s[keep]^2; k = length(v)
  Vt = sum(v); sp2 = sum(v * s2) / Vt
  C  = 1 + (sum(1 / v) - 1 / Vt) / (3 * (k - 1))
  stat = (Vt * log(sp2) - sum(v * log(s2))) / C
  data.frame(chisq = stat, df = k - 1,
             p = pchisq(stat, k - 1, lower.tail = FALSE),
             sd_ratio_max_min = sqrt(max(s2) / min(s2)))
}

## -----------------------------------------------------------------------------
## Comparators
## -----------------------------------------------------------------------------

## Naive weighted regression of stratum means (what lm(mean ~ ..., weights = n)
## returns): same point estimates, but sigma^2 is estimated from SSE1 alone with
## M - K df. Not estimable when the model is saturated (M = K).
naive_wls <- function(formula, data, n = "n", mean = "mean", sd = "sd") {
  p   = .slm_prepare(formula, data, n, mean, sd)
  X   = p$X; ni = p$ni; ybar = p$ybar
  M   = nrow(X); K = ncol(X)
  V   = solve(crossprod(X, X * ni))
  beta = drop(V %*% crossprod(X, ni * ybar))
  df  = M - K
  if (df < 1) {
    se = rep(NA_real_, K)
  } else {
    sigma2 = sum(ni * drop(ybar - X %*% beta)^2) / df
    se = sqrt(diag(V) * sigma2)
  }
  tval = beta / se
  data.frame(term = colnames(X), estimate = beta, se = se, t = tval,
             p = if (df >= 1) 2 * pt(-abs(tval), df) else NA_real_,
             df = df, row.names = NULL)
}

## Collapse strata over all covariates except `by`, returning group n, mean, SD.
collapse_strata <- function(data, by, n = "n", mean = "mean", sd = "sd") {
  data = data[data[[n]] > 0, , drop = FALSE]
  s = data[[sd]]; s[data[[n]] == 1] = 0
  g = interaction(data[by], drop = TRUE)
  out = do.call(rbind, lapply(split(seq_len(nrow(data)), g), function(idx) {
    ni = data[[n]][idx]; yb = data[[mean]][idx]; si = s[idx]
    Ng = sum(ni); mg = sum(ni * yb) / Ng
    ss = sum((ni - 1) * si^2) + sum(ni * (yb - mg)^2)
    res = data[idx[1], by, drop = FALSE]
    res$n = Ng; res$mean = mg; res$sd = if (Ng > 1) sqrt(ss / (Ng - 1)) else NA
    res
  }))
  rownames(out) = NULL
  out
}

## Unadjusted two-sample comparison from summary data (crude analysis used in
## prior studies). Returns pooled-variance and Welch t-tests for group2 - group1.
crude_ttest <- function(data, group, n = "n", mean = "mean", sd = "sd", level = 0.95) {
  g = collapse_strata(data, group, n, mean, sd)
  if (nrow(g) != 2) stop("`group` must define exactly two groups.")
  n1 = g$n[1]; n2 = g$n[2]; m1 = g$mean[1]; m2 = g$mean[2]
  v1 = g$sd[1]^2; v2 = g$sd[2]^2
  diff = m2 - m1
  # pooled
  df_p = n1 + n2 - 2
  se_p = sqrt(((n1 - 1) * v1 + (n2 - 1) * v2) / df_p * (1 / n1 + 1 / n2))
  # Welch
  se_w = sqrt(v1 / n1 + v2 / n2)
  df_w = se_w^4 / ((v1 / n1)^2 / (n1 - 1) + (v2 / n2)^2 / (n2 - 1))
  mk = function(lab, se, df) {
    tv = diff / se; q = qt(1 - (1 - level) / 2, df)
    data.frame(method = lab, estimate = diff, se = se, t = tv, df = df,
               p = 2 * pt(-abs(tv), df), lower = diff - q * se, upper = diff + q * se)
  }
  out <- rbind(mk("t-test (pooled)", se_p, df_p), mk("t-test (Welch)", se_w, df_w))
  attr(out, "groups") <- g
  out
}
