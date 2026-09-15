
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(statmod)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a
lapply_pb <- if (requireNamespace("pbapply", quietly = TRUE)) pbapply::pblapply else lapply

if(!require(acrt))devtools::install_github("samuel-watson/acrt")
source("xo_model_builder.R")
set.seed(20260902)

# -----------------------------------------------------------------------------
# 0. Settings
# -----------------------------------------------------------------------------
NSIM   <- 1000
ALPHA  <- 0.05
Z_CRIT <- qnorm(1 - ALPHA / 2)

## ---- Outcome ---------------------------------------------------------------
FAMILY     <- "binomial"         # "binomial" or "gaussian"
BASELINE_P <- 0.15               # PLANNED control-arm risk
DELTA      <- if (FAMILY == "binomial") -0.03 else 0.25   # risk difference / mean diff

## ---- Design ----------------------------------------------------------------
K1     <- 15L                    # stage 1 clusters PER SEQUENCE
M11    <- 75L                    # cluster size, period 0
M12    <- 75L                    # cluster size, period 1
T2     <- 2L                     # stage 2 periods: the second full crossover
K2MAX  <- 0L                     # max NEW clusters per sequence
M2GRID <- seq(30L, 80L, by = 10L)
T_MAX  <- 2L + T2

RHO   <- 300                     # cost per cluster
RHO_T <- 30                      # cost per period

## ---- Correlation model -----------------------------------------------------
ICC_P <- 0.10                    # planning ICC (see ICC_SCALE)
CAC_P <- 0.80                    # planning CAC = lag-1 cluster autocorrelation

## "ar1"  -> Corr(b_it, b_it') = cac^|t-t'|   (what gr(cl)*ar0(t) fits)
## "exch" -> Corr(b_it, b_it') = cac          (u_i + v_it, the parallel script)
COR_STRUCT_TRUE    <- "ar1"
COR_STRUCT_ASSUMED <- "ar1"
ICC_SCALE <- "marginal"
CONDITIONAL_TARGET <- TRUE
PERIOD_EFFECTS <- rep(0, T_MAX)  # secular trend on the link/mean scale

## ---- Interim behaviour -----------------------------------------------------
KNOWN_THETA       <- FALSE
THETA_RULE        <- "plugin"    # "plugin", "planned", "conservative"
ESTIMATE_CAC      <- FALSE
ESTIMATE_BASELINE <- TRUE                                ### NEW ###
ESTIMATE_ICC      <- TRUE
THETA_HAT_FIELDS  <- c("icc", "cac", "baseline")   ### CHECK ### does
ICC_FLOOR <- 1e-4
CAC_FLOOR <- 0
CAC_CEIL  <- 0.99
P0_FLOOR  <- 0.01
A_FLOOR   <- 1e-6
SMALL_SAMPLE <- FALSE
PROJECTION   <- "full"           # "full" or "stagewise"
DF2_RULE     <- "difference"     # "difference" or "package"
METHOD       <- "lambda"

## Tolerance on |z1(full decomposition) - z1(stage 1 only)|. Exactly zero for
## gaussian; about 1e-4 for binomial, because the GLM working weight
## 1/(m p(1-p)) differs between treated and control cluster-periods and so
## breaks the sequence-swap symmetry of Section 1b by a small amount.
Z1_GAP_TOL <- if (FAMILY == "binomial") 1e-2 else 1e-7          ### CHANGED ###

## ---- Scenarios -------------------------------------------------------------
SCENARIOS <- CJ(hypothesis = c("H0", "H1"),
                icc_true   = c(0.05, 0.10, 0.15),
                p0_true    = c(0.12, 0.15, 0.18), sorted = FALSE)
SCENARIOS[, delta_true := ifelse(hypothesis == "H1", DELTA, 0)]
SCENARIOS[, cac_true := CAC_P]

# -----------------------------------------------------------------------------
# 1. Score decomposition
#
# The manuscript defines the projection stage-wise,
#   x_tilde_t = x_t - X_t (X_t' Sigma_tt^{-1} X_t)^{-1} X_t' Sigma_tt^{-1} x_t,
# whereas this projects once against the full Sigma. Section 1a implements the
# stage-wise version; Section 1b says when they agree.
# -----------------------------------------------------------------------------
efficient_score_decomposition <- function(X, V, idx1, idx2, j = 2, r = NULL) {
  X_trt  <- X[, j, drop = FALSE]
  X_nuis <- X[, -j, drop = FALSE]
  R <- chol(V)
  V_inv_X_nuis <- backsolve(R, backsolve(R, X_nuis, transpose = TRUE))
  V_inv_X_trt  <- backsolve(R, backsolve(R, X_trt,  transpose = TRUE))
  I_nuis_inv <- solve(crossprod(X_nuis, V_inv_X_nuis))
  cross_info <- crossprod(X_nuis, V_inv_X_trt)
  X_trt_tilde <- X_trt - X_nuis %*% (I_nuis_inv %*% cross_info)
  
  X1 <- X_trt_tilde[idx1, , drop = FALSE]
  X2 <- X_trt_tilde[idx2, , drop = FALSE]
  V11 <- V[idx1, idx1, drop = FALSE]
  V12 <- V[idx1, idx2, drop = FALSE]
  V22 <- V[idx2, idx2, drop = FALSE]
  R11 <- chol(V11)
  V11_inv_X1 <- backsolve(R11, backsolve(R11, X1, transpose = TRUE))
  I1_eff <- as.numeric(crossprod(X1, V11_inv_X1))
  V11_inv_V12 <- backsolve(R11, backsolve(R11, V12, transpose = TRUE))
  S   <- V22 - crossprod(V12, V11_inv_V12)
  V21 <- t(V12)
  X2_cond <- X2 - V21 %*% V11_inv_X1
  R_S <- chol(S)
  S_inv_X2 <- backsolve(R_S, backsolve(R_S, X2_cond, transpose = TRUE))
  I2_eff <- as.numeric(crossprod(X2_cond, S_inv_X2))
  
  out <- list(I1_eff = I1_eff, I2_eff = I2_eff, I_eff = I1_eff + I2_eff,
              X_trt_tilde = X_trt_tilde)
  
  if (!is.null(r)) {
    r1 <- r[idx1]; r2 <- r[idx2]
    V11_inv_r1 <- backsolve(R11, backsolve(R11, r1, transpose = TRUE))
    U1 <- as.numeric(crossprod(X1, V11_inv_r1))
    r2_cond  <- r2 - V21 %*% V11_inv_r1
    S_inv_r2 <- backsolve(R_S, backsolve(R_S, r2_cond, transpose = TRUE))
    U2_cond  <- as.numeric(crossprod(X2_cond, S_inv_r2))
    out$U1 <- U1; out$U2_cond <- U2_cond
    out$z1 <- U1 / sqrt(I1_eff)
    out$z2 <- U2_cond / sqrt(I2_eff)
  }
  out
}

.GH <- new.env(parent = emptyenv())

gh_nodes <- function(nodes = 32L) {
  key <- as.character(nodes)
  if (is.null(.GH[[key]])) {
    g <- statmod::gauss.quad(nodes, kind = "hermite")
    .GH[[key]] <- list(x = g$nodes, w = g$weights / sqrt(pi))
  }
  .GH[[key]]
}

## E_b[expit(eta + b)], b ~ N(0, A). Vectorised over eta.
## A not finite / <= 0 => the conditional mean, i.e. the old behaviour.
marginal_p <- function(eta, A, nodes = 32L) {
  if (!is.finite(A) || A <= 0) return(plogis(eta))
  g <- gh_nodes(nodes)
  P <- plogis(outer(as.numeric(eta), g$x * sqrt(2 * A), "+"))
  as.numeric(P %*% g$w)
}

## d/d eta of the above = E_b[p(1-p)]. Reduces to p(1-p) when A is absent.
marginal_dp <- function(eta, A, nodes = 32L) {
  if (!is.finite(A) || A <= 0) { p <- plogis(eta); return(p * (1 - p)) }
  g <- gh_nodes(nodes)
  P <- plogis(outer(as.numeric(eta), g$x * sqrt(2 * A), "+"))
  as.numeric((P * (1 - P)) %*% g$w)
}

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
efficient_score_decomposition_sw <- function(X, V, idx1, idx2, j = 2, r = NULL,
                                             stage2_sigma = c("marginal", "conditional")) {
  stage2_sigma <- match.arg(stage2_sigma)
  
  project <- function(idx, Sig) {
    Xt <- X[idx, , drop = FALSE]
    xt <- Xt[, j, drop = FALSE]
    Xn <- Xt[, -j, drop = FALSE]
    keep <- apply(Xn, 2, function(cc) any(abs(cc) > 0))
    Xn <- Xn[, keep, drop = FALSE]
    if (ncol(Xn) == 0L) return(xt)
    Rc <- chol(Sig)
    SiXn <- backsolve(Rc, backsolve(Rc, Xn, transpose = TRUE))
    Sixt <- backsolve(Rc, backsolve(Rc, xt, transpose = TRUE))
    xt - Xn %*% solve(crossprod(Xn, SiXn), crossprod(Xn, Sixt))
  }
  
  V11 <- V[idx1, idx1, drop = FALSE]
  V12 <- V[idx1, idx2, drop = FALSE]
  V22 <- V[idx2, idx2, drop = FALSE]
  R11 <- chol(V11)
  S   <- V22 - crossprod(V12, backsolve(R11, backsolve(R11, V12, transpose = TRUE)))
  
  X1 <- project(idx1, V11)
  X2 <- project(idx2, if (stage2_sigma == "marginal") V22 else S)
  
  V11_inv_X1 <- backsolve(R11, backsolve(R11, X1, transpose = TRUE))
  I1_eff <- as.numeric(crossprod(X1, V11_inv_X1))
  X2_cond <- X2 - t(V12) %*% V11_inv_X1
  R_S <- chol(S)
  I2_eff <- as.numeric(crossprod(X2_cond,
                                 backsolve(R_S, backsolve(R_S, X2_cond, transpose = TRUE))))
  
  out <- list(I1_eff = I1_eff, I2_eff = I2_eff, I_eff = I1_eff + I2_eff)
  if (!is.null(r)) {
    r1 <- r[idx1]; r2 <- r[idx2]
    V11_inv_r1 <- backsolve(R11, backsolve(R11, r1, transpose = TRUE))
    out$U1 <- as.numeric(crossprod(X1, V11_inv_r1))
    r2_cond <- r2 - t(V12) %*% V11_inv_r1
    out$U2_cond <- as.numeric(crossprod(X2_cond,
                                        backsolve(R_S, backsolve(R_S, r2_cond, transpose = TRUE))))
    out$z1 <- out$U1 / sqrt(I1_eff)
    out$z2 <- out$U2_cond / sqrt(I2_eff)
  }
  out
}

score_decomp <- function(X, V, idx1, idx2, j, r = NULL) {
  if (PROJECTION == "stagewise")
    efficient_score_decomposition_sw(X, V, idx1, idx2, j = j, r = r)
  else
    efficient_score_decomposition(X, V, idx1, idx2, j = j, r = r)
}


marginal_score_z <- function(y, X, V, j = 2) {
  Xt <- X[, j, drop = FALSE]; Xn <- X[, -j, drop = FALSE]
  R  <- chol(V)
  ViXn <- backsolve(R, backsolve(R, Xn, transpose = TRUE))
  ViXt <- backsolve(R, backsolve(R, Xt, transpose = TRUE))
  Xtil <- Xt - Xn %*% solve(crossprod(Xn, ViXn), crossprod(Xn, ViXt))
  Viy  <- backsolve(R, backsolve(R, y, transpose = TRUE))
  bn   <- solve(crossprod(Xn, ViXn), crossprod(Xn, Viy))
  r    <- as.numeric(y - Xn %*% bn)
  Vir  <- backsolve(R, backsolve(R, r, transpose = TRUE))
  ViXtil <- backsolve(R, backsolve(R, Xtil, transpose = TRUE))
  as.numeric(crossprod(Xtil, Vir)) / sqrt(as.numeric(crossprod(Xtil, ViXtil)))
}

null_residual <- function(y, X, V, j = 2) {
  Xn <- X[, -j, drop = FALSE]
  R  <- chol(V)
  ViXn <- backsolve(R, backsolve(R, Xn, transpose = TRUE))
  Viy  <- backsolve(R, backsolve(R, y,  transpose = TRUE))
  as.numeric(y - Xn %*% solve(crossprod(Xn, ViXn), crossprod(Xn, Viy)))
}

t_to_z <- function(tstat, df) {
  if (!SMALL_SAMPLE) return(tstat)
  qnorm(pt(tstat, df = df))
}

## Degrees of freedom: cluster-periods minus fixed effects.
n_periods_fun <- function(t2) 2L + t2
df1_fun <- function() 4L * K1 - 3L
df_full_fun <- function(k2, t2) {
  ncp <- 4L * K1 + t2 * 2L * K1 + t2 * 2L * k2
  ncp - (1L + (n_periods_fun(t2) - 1L) + 1L)
}
## The builder returns df_s2 = t2*2*(k1+k2) - t2 - 1, one LESS than
## df_full - df_s1. Both offered; they cannot both be right.      ### CHECK ###
df2_fun <- function(k2, t2) {
  if (DF2_RULE == "package") t2 * 2L * (K1 + k2) - t2 - 1L
  else df_full_fun(k2, t2) - df1_fun()
}


# -----------------------------------------------------------------------------
# 2. Design layout, scales, and the data generating model
# -----------------------------------------------------------------------------

## Treatment indicator, reproducing crossover_model_builder() exactly:
##   stage 1  sq==0 treated at t==1;  sq==1 treated at t==0
##   stage 2  sq==1 treated at t==2;  sq==0 treated at t==3
trt_rule <- function(sq, t) as.integer(xor(sq == 0L, t %% 2L == 0L))

## Periods observed by clusters recruited AT THE INTERIM. The builder always
## gives them periods 2 and 3 even when t2 == 1, while n_s2 charges them for
## one. "package" reproduces that; it does not bite at t2 = 2.
NEW_CLUSTER_PERIODS <- "package"
new_periods <- function(t2) {
  if (NEW_CLUSTER_PERIODS == "package") 2L:3L else 2L:(1L + t2)
}

design_cp <- function(k1, k2, t2, m11, m12, m2) {
  cl1  <- seq_len(2L * k1)
  seq1 <- rep(0:1, each = k1)
  s1 <- rbind(
    data.table(cl = cl1, sq = seq1, t = 0L, m = m11, stage = 1L),
    data.table(cl = cl1, sq = seq1, t = 1L, m = m12, stage = 1L))
  s2 <- rbindlist(lapply(2L:(1L + t2), function(u)
    data.table(cl = cl1, sq = seq1, t = u, m = m2, stage = 2L)))
  d <- rbind(s1, s2)
  if (k2 > 0L) {
    cl2  <- 2L * k1 + seq_len(2L * k2)
    seq2 <- rep(0:1, each = k2)
    s2n <- rbindlist(lapply(new_periods(t2), function(u)
      data.table(cl = cl2, sq = seq2, t = u, m = m2, stage = 2L)))
    d <- rbind(d, s2n)
  }
  d[, trt := trt_rule(sq, t)]
  d[]
}

## Design matrix, matching ~ factor(t) + trt with treatment LAST.
make_X <- function(cpm) {
  tt <- sort(unique(cpm$t))
  D <- vapply(tt[-1], function(u) as.numeric(cpm$t == u), numeric(nrow(cpm)))
  X <- cbind(1, D, cpm$trt)
  colnames(X) <- c("(Intercept)", paste0("t", tt[-1]), "trt")
  X
}

cor_periods <- function(lag, cac, struct) {
  if (struct == "ar1") cac^lag else ifelse(lag == 0, 1, cac)
}

## -----------------------------------------------------------------------------
## 2a. Scales
##
## gaussian: Var(y_ijt) = 1, cluster variance = icc, individual = 1 - icc.
## binomial: everything on the LINK scale, matching the package, which builds
##   V = W^{-1} + Z D Z' with W the GLM working weights. So the cluster-period
##   diagonal is 1/(m p(1-p)), NOT p(1-p)/m, and the cluster variance is
##   A = icc/(1-icc)/(p0(1-p0)) under ICC_SCALE = "marginal".     ### CHANGED ###
## -----------------------------------------------------------------------------
A_PLAN <- if (FAMILY == "binomial")
  ICC_P / (1 - ICC_P) / (BASELINE_P * (1 - BASELINE_P)) else NA_real_

A_from_icc <- function(icc, p0) {
  if (ICC_SCALE == "link") A_PLAN * (icc / ICC_P)
  else icc / (1 - icc) / (p0 * (1 - p0))
}
icc_from_A <- function(A, p0) {
  if (ICC_SCALE == "link") ICC_P * (A / A_PLAN)
  else { q <- A * p0 * (1 - p0); q / (1 + q) }
}

## Link-scale effect implied by a risk difference at a given control risk.
## This is the crux of baseline re-estimation: the risk-difference target is
## fixed by the protocol, but the link-scale target it implies MOVES:
##   p0 = 0.12 -> -0.3212,  0.15 -> -0.2578,  0.18 -> -0.2183.      ### NEW ###
b1_link_at <- function(p0, rd = DELTA) qlogis(p0 + rd) - qlogis(p0)

vcomp <- function(icc) list(s2c = icc, s2e = 1 - icc)

build_V <- function(cpm, icc, cac, baseline = BASELINE_P,
                    struct = COR_STRUCT_ASSUMED, family = FAMILY,
                    b1_link = NULL) {
  same_cl <- outer(cpm$cl, cpm$cl, "==")
  lag <- abs(outer(cpm$t, cpm$t, "-"))
  if (family == "binomial") {
    A <- A_from_icc(icc, baseline)
    V <- A * cor_periods(lag, cac, struct) * same_cl
    lk <- calibrate_link(baseline, DELTA, A)
    b1 <- if (is.null(b1_link)) lk$b1 else b1_link      ### CHANGED ###
    p <- plogis(lk$eta0 + PERIOD_EFFECTS[cpm$t + 1L] + b1 * cpm$trt)
    diag(V) <- diag(V) + 1 / (cpm$m * p * (1 - p))
  } else {
    s <- vcomp(icc)
    V <- s$s2c * cor_periods(lag, cac, struct) * same_cl
    diag(V) <- diag(V) + s$s2e / cpm$m
  }
  V
}

analysis_y <- function(cpm, V, X, baseline, A = NA_real_,
                       maxit = 50L, tol = 1e-10) {
  if (FAMILY != "binomial") return(cpm$ybar)
  A_use <- if (isTRUE(CONDITIONAL_TARGET) && is.finite(A) && A > 0) A else NA_real_
  sc  <- if (is.finite(A_use)) sqrt(1 + 0.346 * A_use) else 1
  eta <- sc * qlogis(baseline) + PERIOD_EFFECTS[cpm$t + 1L]
  
  R     <- chol(V)
  ViX   <- backsolve(R, backsolve(R, X, transpose = TRUE))
  XtViX <- crossprod(X, ViX)
  
  conv <- FALSE
  for (i in seq_len(maxit)) {
    mu <- marginal_p(eta, A_use)
    d  <- pmax(marginal_dp(eta, A_use), 1e-10)
    z  <- eta + (cpm$ybar - mu) / d
    Viz <- backsolve(R, backsolve(R, z, transpose = TRUE))
    eta_new <- as.numeric(X %*% solve(XtViX, crossprod(X, Viz)))
    if (max(abs(eta_new - eta)) < tol) { eta <- eta_new; conv <- TRUE; break }
    eta <- eta_new
  }
  mu <- marginal_p(eta, A_use)
  d  <- pmax(marginal_dp(eta, A_use), 1e-10)
  structure(eta + (cpm$ybar - mu) / d, converged = conv)
}

analysis_y_null <- function(cpm, V, X, j, baseline, A = NA_real_,
                            maxit = 50L, tol = 1e-10) {
  if (FAMILY != "binomial") return(cpm$ybar)
  A_use <- if (isTRUE(CONDITIONAL_TARGET) && is.finite(A) && A > 0) A else NA_real_
  Xn  <- X[, -j, drop = FALSE]
  sc  <- if (is.finite(A_use)) sqrt(1 + 0.346 * A_use) else 1
  eta <- sc * qlogis(baseline) + PERIOD_EFFECTS[cpm$t + 1L]
  
  R      <- chol(V)
  ViXn   <- backsolve(R, backsolve(R, Xn, transpose = TRUE))
  XnViXn <- crossprod(Xn, ViXn)
  
  conv <- FALSE
  for (i in seq_len(maxit)) {
    mu <- marginal_p(eta, A_use)
    d  <- pmax(marginal_dp(eta, A_use), 1e-10)
    z  <- eta + (cpm$ybar - mu) / d
    Viz <- backsolve(R, backsolve(R, z, transpose = TRUE))
    eta_new <- as.numeric(Xn %*% solve(XnViXn, crossprod(Xn, Viz)))
    if (max(abs(eta_new - eta)) < tol) { eta <- eta_new; conv <- TRUE; break }
    eta <- eta_new
  }
  mu <- marginal_p(eta, A_use)
  d  <- pmax(marginal_dp(eta, A_use), 1e-10)
  structure(eta + (cpm$ybar - mu) / d, converged = conv)
}

b1_hat_stage1 <- function(delta_true = DELTA, icc_true = ICC_P,
                          cac_true = CAC_P, p0_true = BASELINE_P,
                          nsim = 500L, seed = 20260902L) {
  set.seed(seed)
  lp  <- link_pars(p0_true, delta_true, icc_true)
  out <- vapply(seq_len(nsim), function(i) {
    eff <- draw_cluster_effects(K1, icc_true, cac_true, p0_true)
    cl1 <- unlist(lapply(0:1, function(s) which(eff$sqv == s)[seq_len(K1)]))
    cpm1 <- sim_periods_cp(eff, cl1, periods = 0:1,
                           m_by_period = c(M11, M12),
                           delta = delta_true, icc = icc_true, p0 = p0_true)
    ## Planning values, so V is deterministic and only sampling varies.
    V1 <- build_V(cpm1, ICC_P, CAC_P, baseline = BASELINE_P, b1_link = lp$b1)
    X1 <- make_X(cpm1)
    j1 <- which(colnames(X1) == "trt")
    yw <- analysis_y(cpm1, V1, X1, BASELINE_P, A = A_from_icc(ICC_P, BASELINE_P))
    if (!isTRUE(attr(yw, "converged"))) return(NA_real_)
    R   <- chol(V1)
    ViX <- backsolve(R, backsolve(R, X1, transpose = TRUE))
    Viy <- backsolve(R, backsolve(R, as.numeric(yw), transpose = TRUE))
    solve(crossprod(X1, ViX), crossprod(X1, Viy))[j1]
  }, numeric(1))
  
  lk <- calibrate_link(BASELINE_P, DELTA, A_PLAN)
  c(mean       = mean(out, na.rm = TRUE),
    mcse       = sd(out, na.rm = TRUE) / sqrt(sum(!is.na(out))),
    n_nonconv  = sum(is.na(out)),
    target_cond = lk$b1,
    target_marg = b1_link_at(BASELINE_P, DELTA),
    ratio      = mean(out, na.rm = TRUE) / lk$b1)
}

## -----------------------------------------------------------------------------
## 2b. Marginal calibration of the binomial DGM                     ### NEW ###
##
## With a cluster random effect of variance A on the logit scale, expit of the
## conditional intercept is NOT the marginal prevalence: E[expit(eta + b)] is
## pulled toward 1/2. At A = 0.87 that is a couple of percentage points, enough
## to bias the interim baseline estimator and make the believed-vs-actual
## conditional power gap look like a bug when it is a scale artefact.
##
## MARGINAL_TARGET = TRUE solves for the conditional (eta0, b1) delivering the
## requested MARGINAL control risk and risk difference. FALSE uses the naive
## conditional values, which is what b1_link_at() returns. Whichever you pick,
## it must agree with convert_to_glm_scale() or the power comparison carries
## the difference. Section 3 prints all three for comparison.      ### CHECK ###
## -----------------------------------------------------------------------------
MARGINAL_TARGET <- TRUE

marginal_p <- function(eta, A, nodes = 32L) {
  gh <- statmod::gauss.quad(nodes, kind = "hermite")
  x <- gh$nodes * sqrt(2 * A); w <- gh$weights / sqrt(pi)
  sum(w * plogis(eta + x))
}

calibrate_link <- function(p0, rd, A) {
  if (!MARGINAL_TARGET || !is.finite(A) || A <= 0)
    return(list(eta0 = qlogis(p0), b1 = if (rd == 0) 0 else b1_link_at(p0, rd)))
  e0 <- uniroot(function(e) marginal_p(e, A) - p0,
                qlogis(p0) + c(-4, 4))$root
  if (rd == 0) return(list(eta0 = e0, b1 = 0))
  e1 <- uniroot(function(e) marginal_p(e, A) - (p0 + rd),
                qlogis(p0 + rd) + c(-4, 4))$root
  list(eta0 = e0, b1 = e1 - e0)
}

## -----------------------------------------------------------------------------
## 2c. Simulation of cluster-period sufficient statistics
##
## The score depends on the data only through cluster-period means, and the
## interim variance-component estimator only through the pooled within-
## cluster-period sum of squares. Both are simulated directly. Individual-level
## data is generated only in the verification section.
## -----------------------------------------------------------------------------
draw_cluster_effects <- function(k_max_per_seq, icc, cac, p0,
                                 struct = COR_STRUCT_TRUE, nT = T_MAX) {
  s2c <- if (FAMILY == "binomial") A_from_icc(icc, p0) else vcomp(icc)$s2c
  R <- s2c * cor_periods(abs(outer(seq_len(nT), seq_len(nT), "-")), cac, struct)
  L <- chol(R + diag(1e-12, nT))
  n_cl <- 2L * k_max_per_seq
  list(B = matrix(rnorm(n_cl * nT), n_cl, nT) %*% L,
       sqv = rep(0:1, each = k_max_per_seq),
       s2c = s2c)
}

sim_periods_cp <- function(eff, cl_use, periods, m_by_period, delta, icc, p0) {
  lk  <- if (FAMILY == "binomial") calibrate_link(p0, delta, eff$s2c) else NULL
  s2e <- if (FAMILY == "binomial") NA_real_ else vcomp(icc)$s2e
  rbindlist(lapply(seq_along(periods), function(a) {
    tt <- periods[a]; m <- m_by_period[a]
    sq <- eff$sqv[cl_use]
    tr <- trt_rule(sq, tt)
    n  <- length(cl_use)
    if (FAMILY == "binomial") {
      eta  <- lk$eta0 + PERIOD_EFFECTS[tt + 1L] + lk$b1 * tr + eff$B[cl_use, tt + 1L]
      ybar <- rbinom(n, m, plogis(eta)) / m
      ss   <- m * ybar * (1 - ybar)
    } else {
      mu   <- PERIOD_EFFECTS[tt + 1L] + delta * tr + eff$B[cl_use, tt + 1L]
      ybar <- mu + rnorm(n, 0, sqrt(s2e / m))
      ss   <- s2e * rchisq(n, df = m - 1L)
    }
    data.table(cl = cl_use, sq = sq, t = tt, trt = tr, m = m,
               ybar = ybar, ss = ss, dfw = m - 1L)
  }))
}

## -----------------------------------------------------------------------------
## 2d. Interim estimation of (icc, cac, baseline) from the first crossover
##
## Centre the cluster-period statistics within each sequence-by-period cell.
## That absorbs the intercept, the period effect, the treatment effect AND any
## sequence effect, at exactly one df per cell. With k clusters per sequence:
##   s2u_hat  = sum_i r_i0 r_i1 / (2(k-1))     lag-1 cluster covariance
##   sigma2_t = sum_i r_it^2    / (2(k-1))     cluster-period variance
## and the within-cluster-period component is subtracted off: estimated from
## the pooled sums of squares for gaussian, and KNOWN given m and p-hat for
## binomial, where it is the working variance 1/(m p(1-p)).
##
## Under AR1 this estimates the lag-1 autocorrelation, exactly the ar0()
## parameter; under exchangeability the common correlation. The two are
## indistinguishable from two periods, which is the point of Section 11.
##
## The baseline comes from the control cluster-periods -- available here and
## not in the parallel design's single stage 1 period.              ### NEW ###
## -----------------------------------------------------------------------------
theta_stage1 <- function(cpm1) {
  k <- cpm1[t == 0L & sq == 0L, .N]
  denom <- 2 * (k - 1)
  
  ctl <- cpm1[trt == 0L]
  p0_hat <- min(max(sum(ctl$ybar * ctl$m) / sum(ctl$m), P0_FLOOR), 1 - P0_FLOOR)
  
  if (FAMILY == "binomial") {
    ## empirical logits with a 0.5 continuity correction, and their KNOWN
    ## within-cluster-period working variance
    d <- copy(cpm1)
    d[, yc := pmin(pmax(ybar, 0.5 / m), 1 - 0.5 / m)]
    d[, `:=`(eta = qlogis(yc), wvar = 1 / (m * yc * (1 - yc)))]
    d[, res := eta - mean(eta), by = .(sq, t)]
    w <- dcast(d, cl + sq ~ t, value.var = "res")
    setnames(w, c("0", "1"), c("r0", "r1"))
    s2u_raw <- sum(w$r0 * w$r1) / denom
    v0 <- sum(w$r0^2) / denom; v1 <- sum(w$r1^2) / denom
    A_raw <- mean(c(v0 - d[t == 0L, mean(wvar)], v1 - d[t == 1L, mean(wvar)]))
    A_hat <- max(A_raw, A_FLOOR)
    s2u   <- max(min(s2u_raw, A_hat), 0)
    icc   <- max(icc_from_A(A_hat, p0_hat), ICC_FLOOR)
    cac   <- if (A_hat <= A_FLOOR) CAC_P else min(max(s2u / A_hat, CAC_FLOOR), CAC_CEIL)
    icc_raw <- icc_from_A(max(A_raw, 1e-12), p0_hat)
    cac_raw <- if (A_raw > 0) s2u_raw / A_raw else NA_real_
  } else {
    s2e <- cpm1[, sum(ss) / sum(dfw)]
    r <- copy(cpm1)[, res := ybar - mean(ybar), by = .(sq, t)]
    w <- dcast(r, cl + sq ~ t, value.var = "res")
    setnames(w, c("0", "1"), c("r0", "r1"))
    s2u_raw <- sum(w$r0 * w$r1) / denom
    v0 <- sum(w$r0^2) / denom; v1 <- sum(w$r1^2) / denom
    m0 <- cpm1[t == 0L, m[1]]; m1 <- cpm1[t == 1L, m[1]]
    s2c_raw <- mean(c(v0 - s2e / m0, v1 - s2e / m1))
    s2u <- max(s2u_raw, 0); s2v <- max(s2c_raw - s2u_raw, 0)
    s2c <- s2u + s2v
    icc <- max(s2c / (s2c + s2e), ICC_FLOOR)
    cac <- if (s2c <= 0) CAC_P else min(max(s2u / s2c, CAC_FLOOR), CAC_CEIL)
    icc_raw <- s2c_raw / (s2c_raw + s2e)
    cac_raw <- if (s2c_raw > 0) s2u_raw / s2c_raw else NA_real_
  }
  
  list(icc = if (ESTIMATE_ICC) icc else ICC_P,
       cac = if (ESTIMATE_CAC) cac else CAC_P,
       baseline = if (ESTIMATE_BASELINE) p0_hat else BASELINE_P,
       icc_raw = icc_raw, cac_raw = cac_raw, p0_raw = p0_hat)
}

apply_theta_rule <- function(est) {
  if (KNOWN_THETA)
    return(list(icc = ICC_P, cac = CAC_P, baseline = BASELINE_P))
  switch(THETA_RULE,
         plugin  = list(icc = est$icc, cac = est$cac, baseline = est$baseline),
         planned = list(icc = ICC_P, cac = CAC_P, baseline = BASELINE_P),
         ## Conservative = adapt only in the direction that makes stage 2 look
         ## LESS informative: icc up, cac DOWN (a crossover profits from a
         ## large cac, the opposite sign to the parallel case), baseline UP (a
         ## higher control risk means a SMALLER link-scale effect for a fixed
         ## risk difference).                                     ### CHECK ###
         conservative = list(icc = max(est$icc, ICC_P),
                             cac = min(est$cac, CAC_P),
                             baseline = max(est$baseline, BASELINE_P)),
         stop("unknown THETA_RULE"))
}


# -----------------------------------------------------------------------------
# 3. Design, weights, boundaries
# -----------------------------------------------------------------------------
design_single <- crossover_crt(
  icc      = ICC_P,
  cac      = CAC_P,
  delta    = DELTA,
  baseline = if (FAMILY == "binomial") BASELINE_P else NULL,
  k1       = K1,
  m11      = M11,
  m12      = M12,
  k2       = 0:K2MAX,
  m2       = M2GRID,
  t2       = T2,
  rho      = RHO,
  rho_t    = RHO_T,
  family   = FAMILY
)


## Now correct rather than dangerous: with the covariance parameterisation
## fixed in both branches of the builder, caching the glmmrBase model is a
## pure speed-up. Against an unpatched build it applies the update branch to
## every model after the first, which the assertion below catches.
design_single$spec$model_builder <- make_crossover_model_fn()

results_single <- adaptive_analysis(design_single, target_power = 0.8,
                                    tol = 0.005, method = METHOD)

## I1_eff is a stage 1 quantity: exactly invariant to the stage 2 design for
## gaussian, invariant to 0.015% for binomial (Section 1b). Anything larger
## means the package build is missing the covariance fix and every number
## below is meaningless.                                         ### CHANGED ###
.I1 <- vapply(results_single$raw$results$models$list, function(m) m$I1_eff, numeric(1))
if (diff(range(.I1)) > 1e-3 * mean(.I1))
  stop(sprintf(paste("I1_eff varies across the stage 2 grid (%.6g to %.6g).",
                     "The builder is passing glm_params$cov_pars straight to",
                     "gr(cl)*ar0(t) somewhere -- see header item (a)."),
               min(.I1), max(.I1)))

.params <- results_single$raw$results$params
W1 <- as.numeric(.params$w1_ref)[1]
W2 <- as.numeric(.params$w2_ref %||% .params$w2 %||% sqrt(1 - W1^2))[1]
if (!is.finite(W2)) W2 <- sqrt(1 - W1^2)
stopifnot(is.finite(W1), is.finite(W2), abs(W1^2 + W2^2 - 1) < 1e-6)

## Boundary on the Z scale. Under H0 the transformed statistic is exactly
## standard normal, so this is qnorm-based and df_s1 affects power only.
C_EFF <- as.numeric(.params$efficacy_boundary %||% (Z_CRIT / W1))[1]
C2    <- as.numeric(.params$c2 %||% Z_CRIT)[1]
if (is.null(.params$c2))
  warning("c2 not found in params; (c1, c2) are solved jointly, so the ",
          "z_{alpha/2} fallback is NOT the calibrated boundary.")

results_binding <- tryCatch(
  adaptive_analysis(design_single, target_power = 0.8, tol = 0.005,
                    method = METHOD, futility = "binding"),
  error = function(e) { warning("binding analysis failed: ", conditionMessage(e)); NULL })

CTX <- list(nonbinding = list(fut = "nonbinding", res = results_single, c2 = C2))
if (!is.null(results_binding)) {
  .pb <- results_binding$raw$results$params
  stopifnot(isTRUE(all.equal(as.numeric(.pb$w1_ref)[1], W1)))
  CTX$binding <- list(fut = "binding", res = results_binding,
                      c2 = as.numeric(.pb$c2 %||% C2)[1])
} else {
  cat("binding futility context skipped; running nonbinding only\n")
}

rules <- get_decision_rules(results_single)
EFF_SIGN <- as.numeric(.params$eff_sign %||% 1)[1]
.f <- rules$z1[!rules$continue & abs(rules$z1) < C_EFF]
CONT <- if (!length(.f)) c(-C_EFF, C_EFF) else
  if (EFF_SIGN > 0) c(max(.f), C_EFF) else c(-C_EFF, min(.f))

cat(sprintf("\nw1 = %.4f, efficacy |z1| > %.3f, c2 = %.4f, df1 = %d\n",
            W1, C_EFF, C2, df1_fun()))
if (FAMILY == "binomial") {
  .lk <- calibrate_link(BASELINE_P, DELTA, A_PLAN)
  cat(sprintf("link scale: naive b1 = %.4f | marginally calibrated b1 = %.4f | package b1 = %s\n",
              b1_link_at(BASELINE_P), .lk$b1,
              format(.params$b1 %||% NA_real_, digits = 4)))
  cat("These must agree, or the power comparison carries the difference.\n")
}

exact_score_stage1 <- function(cpm1, psi, A, cac, nodes = 24L) {
  
  d <- data.table::copy(cpm1)
  data.table::setorder(d, cl, t)
  ncl <- data.table::uniqueN(d$cl)
  stopifnot(nrow(d) == 2L * ncl, all(sort(unique(d$t)) == c(0L, 1L)))
  
  Y  <- matrix(round(d$ybar * d$m), ncl, 2L, byrow = TRUE)
  M  <- matrix(d$m,                 ncl, 2L, byrow = TRUE)
  TR <- matrix(d$trt,               ncl, 2L, byrow = TRUE)
  
  ## 2-D Gauss-Hermite on b ~ N(0, Sb):  b = sqrt(2) C x,  C C' = Sb
  g  <- statmod::gauss.quad(nodes, kind = "hermite")
  Sb <- A * matrix(c(1, cac, cac, 1), 2L, 2L)
  C  <- t(chol(Sb))
  
  gr <- as.matrix(expand.grid(j = seq_len(nodes), k = seq_len(nodes)))
  B  <- C %*% (rbind(g$nodes[gr[, "j"]], g$nodes[gr[, "k"]]) * sqrt(2))
  W  <- g$weights[gr[, "j"]] * g$weights[gr[, "k"]] / pi
  
  p0 <- plogis(psi[1] + B[1, ])
  p1 <- plogis(psi[2] + B[2, ])
  
  ## log f(y_i | b) at every node: ncl x nq
  LL <- outer(Y[, 1], log(p0)) + outer(M[, 1] - Y[, 1], log1p(-p0)) +
    outer(Y[, 2], log(p1)) + outer(M[, 2] - Y[, 2], log1p(-p1))
  
  Wt  <- sweep(exp(LL - apply(LL, 1, max)), 2, W, "*")
  Lk  <- rowSums(Wt)
  Ep0 <- as.numeric(Wt %*% p0) / Lk          # E[p_i0 | y_i]
  Ep1 <- as.numeric(Wt %*% p1) / Lk
  
  ## x_tilde = trt - 1/2; in a crossover trt_i0 = 1 - trt_i1, so this is a
  ## within-cluster contrast of the two periods' Pearson-type residuals.
  Ui <- (TR[, 1] - 0.5) * (Y[, 1] - M[, 1] * Ep0) +
    (TR[, 2] - 0.5) * (Y[, 2] - M[, 2] * Ep1)
  
  list(Ui = Ui,
       U  = sum(Ui),
       I_op = sum(Ui^2),                      # outer-product information
       z  = sum(Ui) / sqrt(sum(Ui^2)))
}


# -----------------------------------------------------------------------------
# 5. One replicate
# -----------------------------------------------------------------------------
## Cost = participants + rho * clusters + rho_t * periods, matching
## crt_design_spec(weights = c(n = 1, clusters = "rho", t = "rho_t")).
## Periods are counted once, not per arm.
trial_cost <- function(n_ind, n_clus, n_per) n_ind + RHO * n_clus + RHO_T * n_per
n_stage1 <- function() 2L * K1 * (M11 + M12)
n_stage2 <- function(k2, m2, t2) {
  2L * K1 * m2 * t2 + 2L * k2 * m2 * length(new_periods(t2))
}

## Non-centrality scale: the response scale for gaussian, the link scale for
## binomial -- AT THE BASELINE IN PLAY. The believed conditional power uses the
## interim estimate; cp_actual uses the truth.                      ### NEW ###
delta_scale <- function(d, p0) {
  if (FAMILY != "binomial") return(d)
  if (d == 0) 0 else b1_link_at(p0, d)
}

link_pars <- function(p0, rd, icc) {
  A <- A_from_icc(icc, p0)
  lk <- calibrate_link(p0, rd, A)
  
  list(
    A = A,
    eta0 = lk$eta0,
    b1 = lk$b1
  )
}

run_once <- function(delta_true, icc_true, cac_true, p0_true, ctx, fut) {
  stopifnot(fut %in% c("binding", "nonbinding"))
  .row <- list(
    decision     = NA_character_,
    reject       = NA,
    z1           = NA_real_,
    z2           = NA_real_,
    z2_trueV     = NA_real_,
    cstat        = NA_real_,
    icc_hat      = NA_real_,
    cac_hat      = NA_real_,
    p0_hat       = NA_real_,
    b1_hat       = NA_real_,
    k2           = NA_integer_,
    m2           = NA_integer_,
    t2           = NA_integer_,
    N            = NA_integer_,
    n_clus       = NA_integer_,
    n_per        = NA_integer_,
    cost         = NA_real_,
    cp           = NA_real_,
    cp_own       = NA_real_,
    cp_pkg_gap   = NA_real_,
    cp_actual    = NA_real_,
    cp_gap       = NA_real_,
    I2_true      = NA_real_,
    I2_ratio     = NA_real_,
    z1_gap       = NA_real_,
    recalibrated = NA,
    z_pooled     = NA_real_,
    inconsistent = NA,
    cf_realised  = NA_real_,
    nonconv      = FALSE
  )
  
  out <- function(...) {
    vals <- list(...)
    extra <- setdiff(names(vals), names(.row))
    if (length(extra))
      stop("run_once: field(s) not in template: ", paste(extra, collapse = ", "))
    bad <- vapply(vals, function(v) is.null(v) || length(v) != 1L, logical(1))
    if (any(bad))
      stop("run_once: field(s) not length 1: ",
           paste(names(vals)[bad], collapse = ", "))
    as.data.table(modifyList(.row, vals))
  }
  
  eff <- draw_cluster_effects(K1 + K2MAX, icc_true, cac_true, p0_true)
  cl1 <- unlist(lapply(0:1, function(s) which(eff$sqv == s)[seq_len(K1)]))
  
  ## ---- stage 1: the first full crossover ----
  cpm1 <- sim_periods_cp(eff, cl1, periods = 0:1,
                         m_by_period = c(M11, M12),
                         delta = delta_true, icc = icc_true, p0 = p0_true)
  
  est   <- theta_stage1(cpm1)
  theta <- apply_theta_rule(est)
  lp_int <- link_pars(theta$baseline, DELTA, theta$icc)
  b1_int <- lp_int$b1
  
  V1 <- build_V(cpm1, theta$icc, theta$cac, baseline = theta$baseline,
                b1_link = b1_int)
  X1 <- make_X(cpm1)
  j1 <- which(colnames(X1) == "trt")
  
  yw1 <- analysis_y(cpm1, V1, X1, theta$baseline, A = lp_int$A)
  if (!isTRUE(attr(yw1, "converged"))) .nonconv <- TRUE
  t1  <- marginal_score_z(yw1, X1, V1, j = j1)
  z1  <- t_to_z(t1, df1_fun())
  R1   <- chol(V1)
  ViX1 <- backsolve(R1, backsolve(R1, X1, transpose = TRUE))
  Viy1 <- backsolve(R1, backsolve(R1, as.numeric(yw1), transpose = TRUE))
  b1_hat_s1 <- solve(crossprod(X1, ViX1), crossprod(X1, Viy1))[j1]
  .nonconv <- FALSE
  dec <- withCallingHandlers(
    interim_analysis(ctx$res, z1_obs = z1,
                     theta_hat = theta[THETA_HAT_FIELDS], verbose = FALSE, override_stop = fut == "nonbinding"),
    warning = function(w) {
      if (grepl("Maximum iterations", conditionMessage(w))) {
        .nonconv <<- TRUE
        invokeRestart("muffleWarning")
      }
    })
  
  stopifnot(isTRUE(all.equal(as.numeric(dec$w1_ref)[1], W1)))
  if (!is.null(dec$c2)) stopifnot(abs(dec$c2 - ctx$c2) < 1e-6)
  
  n1  <- n_stage1()
  nc1 <- 2L * K1
  
  if ((dec$decision != "continue" & fut == "binding") | (dec$decision == "efficacy_stop" & fut == "nonbinding")) {
    return(out(
      decision    = dec$decision,
      reject      = dec$decision == "efficacy_stop",
      z1          = z1,
      icc_hat     = theta$icc,
      cac_hat     = theta$cac,
      p0_hat      = theta$baseline,
      b1_hat      = b1_hat_s1,
      k2 = 0L, m2 = 0L, t2 = 0L,
      N           = n1,
      n_clus      = nc1,
      n_per       = 2L,
      cost        = trial_cost(n1, nc1, 2L),
      cf_realised = if (dec$decision == "futility_stop") z1 else NA_real_,
      nonconv     = .nonconv))
  }
  
  k2 <- as.integer(dec$stage2_design$k2)
  m2 <- as.integer(dec$stage2_design$m2)
  t2 <- as.integer(dec$stage2_design$t2 %||% T2)
  
  ## ---- stage 2: the second full crossover ----
  cl_new <- if (k2 > 0L)
    unlist(lapply(0:1, function(s) setdiff(which(eff$sqv == s), cl1)[seq_len(k2)]))
  else integer(0)
  
  per_cont <- 2L:(1L + t2)
  cpm2 <- sim_periods_cp(eff, cl1, periods = per_cont,
                         m_by_period = rep(m2, length(per_cont)),
                         delta = delta_true, icc = icc_true, p0 = p0_true)
  if (k2 > 0L) {
    per_new <- new_periods(t2)
    cpm2 <- rbind(cpm2,
                  sim_periods_cp(eff, cl_new, periods = per_new,
                                 m_by_period = rep(m2, length(per_new)),
                                 delta = delta_true, icc = icc_true, p0 = p0_true))
  }
  
  cpm  <- rbind(cpm1, cpm2)
  idx1 <- which(cpm$t <= 1L); idx2 <- which(cpm$t >= 2L)
  V  <- build_V(cpm, theta$icc, theta$cac, baseline = theta$baseline,
                b1_link = b1_int)
  X  <- make_X(cpm)
  jt <- which(colnames(X) == "trt")
  r <- null_residual(analysis_y(cpm, V, X, theta$baseline, A = lp_int$A), X, V, j = jt)
  sc <- score_decomp(X, V, idx1, idx2, j = jt, r = r)
  z1_gap <- t_to_z(sc$z1, df1_fun()) - z1
  lp_true <- link_pars(p0_true, delta_true, icc_true)
  b1_true <- lp_true$b1
  V_true  <- build_V(cpm, icc_true, cac_true, baseline = p0_true,
                     struct = COR_STRUCT_TRUE, b1_link = b1_true)
  sc_true <- score_decomp(X, V_true, idx1, idx2, j = jt)
  r_true <- null_residual(analysis_y(cpm, V_true, X, p0_true, A = lp_true$A),
                          X, V_true, j = jt)
  
  sc_true_obs <- score_decomp(
    X, V_true, idx1, idx2,
    j = jt, r = r_true
  )
  
  z2_trueV <- sc_true_obs$z2
  cp_actual <- pnorm((W1 * z1 - ctx$c2) / W2 + b1_true * sqrt(sc_true$I2_eff))
  z2 <- t_to_z(sc$z2, df2_fun(k2, t2))
  I2_ratio <- if (is.null(dec$I2_eff_updated)) NA_real_
  else sc$I2_eff / dec$I2_eff_updated
  cstat  <- W1 * z1 + W2 * z2
  n_tot  <- n1 + n_stage2(k2, m2, t2)
  nc_tot <- 2L * (K1 + k2)
  np_tot <- n_periods_fun(t2)
  z_pooled <- (sqrt(sc$I1_eff) * z1 + sqrt(sc$I2_eff) * z2) /
    sqrt(sc$I1_eff + sc$I2_eff)
  
  cp_believed <- if (is.null(dec$conditional_power)) NA_real_
  else dec$conditional_power
  cp_own <- pnorm((W1 * z1 - ctx$c2) / W2 + b1_int * sqrt(sc$I2_eff))
  
  out(decision     = dec$decision,
      reject       = abs(cstat) > ctx$c2 | dec$decision == "efficacy_stop",
      z1           = z1,
      z2           = z2,
      z2_trueV     = z2_trueV,
      cstat        = cstat,
      icc_hat      = theta$icc,
      cac_hat      = theta$cac,
      p0_hat       = theta$baseline,
      b1_hat       = b1_hat_s1,
      k2 = k2, m2 = m2, t2 = t2,
      N            = n_tot,
      n_clus       = nc_tot,
      n_per        = np_tot,
      cost         = trial_cost(n_tot, nc_tot, np_tot),
      cp           = cp_believed,
      cp_own       = cp_own,
      cp_pkg_gap   = cp_believed - cp_own,
      cp_actual    = cp_actual,
      cp_gap       = cp_believed - cp_actual,
      I2_true      = sc_true$I2_eff,
      I2_ratio     = I2_ratio,
      z1_gap       = z1_gap,
      recalibrated = isTRUE(dec$recalibrated),
      z_pooled     = z_pooled,
      inconsistent = abs(cstat) > ctx$c2 & sign(cstat) != sign(z_pooled),
      nonconv      = .nonconv)
}

## run_scenario() was called in Sections 11-14 of the parallel script but never
## defined there.
# run_scenario <- function(delta_true, icc_true, cac_true = CAC_P,
#                          p0_true = BASELINE_P, ctx = CTX$nonbinding,
#                          nsim = NSIM, fut = "nonbinding") {
#   rbindlist(lapply_pb(seq_len(nsim), function(i)
#     run_once(delta_true, icc_true, cac_true, p0_true, ctx, fut)))
# }

# =============================================================================

## -----------------------------------------------------------------------------
## Posterior cluster-period risks under an AR1 cluster random effect.
##
## b_1 ~ N(0, A);  b_t | b_{t-1} ~ N(cac * b_{t-1}, A(1 - cac^2))
##
## The b_it form a Markov chain, so the T-dimensional integral factorises into T
## one-dimensional quadratures via forward-backward. Cost is linear in T.
##
##   Y, M, TRT, ETA : ncl x T matrices (counts, denominators, treatment,
##                    NULL linear predictor with no treatment term)
##   obs            : ncl x T logical, FALSE where a cluster is unobserved
##                    (new clusters in stage 1 periods)
##
## Returns ncl x T matrix of E[p_it | y_i].
## -----------------------------------------------------------------------------
posterior_p_ar1 <- function(Y, M, ETA, A, cac, obs = NULL,
                            ngrid = 80L, span = 6) {
  ncl <- nrow(Y); Tn <- ncol(Y)
  if (is.null(obs)) obs <- matrix(TRUE, ncl, Tn)
  
  b <- seq(-span * sqrt(A), span * sqrt(A), length.out = ngrid)
  h <- b[2] - b[1]
  
  ## transition K[j,k] = p(b_t = b_j | b_{t-1} = b_k)
  sd_i <- sqrt(max(A * (1 - cac^2), 1e-12))
  K <- outer(b, cac * b, function(u, v) dnorm(u, v, sd_i))
  prior <- dnorm(b, 0, sqrt(A))
  
  ## G[[t]] : ncl x ngrid, the binomial likelihood at each grid point
  G <- lapply(seq_len(Tn), function(t) {
    P <- plogis(outer(ETA[, t], b, "+"))
    L <- matrix(1, ncl, ngrid)
    o <- obs[, t]
    if (any(o))
      L[o, ] <- exp(Y[o, t] * log(P[o, , drop = FALSE]) +
                      (M[o, t] - Y[o, t]) * log1p(-P[o, , drop = FALSE]))
    L
  })
  
  nrm <- function(X) X / pmax(rowSums(X), .Machine$double.xmin)
  
  ## forward
  alpha <- vector("list", Tn)
  alpha[[1]] <- nrm(G[[1]] * rep(prior, each = ncl))
  for (t in 2:Tn)
    alpha[[t]] <- nrm(G[[t]] * (alpha[[t - 1]] %*% t(K)) * h)
  
  ## backward
  beta <- vector("list", Tn)
  beta[[Tn]] <- matrix(1, ncl, ngrid)
  for (t in (Tn - 1):1)
    beta[[t]] <- nrm((G[[t + 1]] * beta[[t + 1]]) %*% K * h)
  
  Ep <- matrix(NA_real_, ncl, Tn)
  for (t in seq_len(Tn)) {
    g <- alpha[[t]] * beta[[t]]
    P <- plogis(outer(ETA[, t], b, "+"))
    Ep[, t] <- rowSums(g * P) / pmax(rowSums(g), .Machine$double.xmin)
  }
  Ep
}

score_contribs <- function(cpm, psi, A, cac, ngrid = 80L) {
  d   <- data.table::copy(cpm)
  cls <- sort(unique(d$cl))
  tt  <- sort(unique(d$t))
  ncl <- length(cls); Tn <- length(tt)
  
  idx <- cbind(match(d$cl, cls), match(d$t, tt))
  Y <- M <- TRT <- matrix(0, ncl, Tn)
  obs <- matrix(FALSE, ncl, Tn)
  Y[idx]   <- round(d$ybar * d$m)
  M[idx]   <- d$m
  TRT[idx] <- d$trt
  obs[idx] <- TRUE
  
  ETA <- matrix(psi[tt + 1L], ncl, Tn, byrow = TRUE)
  Ep  <- posterior_p_ar1(Y, M, ETA, A, cac, obs, ngrid)
  
  Xt <- (TRT - 0.5) * obs                      # x_tilde, zero where unobserved
  setNames(rowSums(Xt * (Y - M * Ep) * obs), as.character(cls))
}


## -----------------------------------------------------------------------------
## Whole-trial replicate using the exact score.

run_once_exact <- function(delta_true, icc_true, cac_true, p0_true, ctx, fut,
                           ngrid = 80L) {
  
  A_true <- A_from_icc(icc_true, p0_true)
  eff <- draw_cluster_effects(K1 + K2MAX, icc_true, cac_true, p0_true)
  cl1 <- unlist(lapply(0:1, function(s) which(eff$sqv == s)[seq_len(K1)]))
  
  cpm1 <- sim_periods_cp(eff, cl1, periods = 0:1, m_by_period = c(M11, M12),
                         delta = delta_true, icc = icc_true, p0 = p0_true)
  
  est   <- theta_stage1(cpm1)
  theta <- apply_theta_rule(est)
  A_hat <- A_from_icc(theta$icc, theta$baseline)
  ## Null linear predictor: no treatment term. A score test evaluates at beta=0
  ## under H1 too. Nuisance fixed from stage 1, as elsewhere.
  psi <- rep(calibrate_link(theta$baseline, 0, A_hat)$eta0, T_MAX) + PERIOD_EFFECTS
  
  U1i <- score_contribs(cpm1, psi, A_hat, theta$cac, ngrid)
  I1  <- sum(U1i^2)
  z1  <- sum(U1i) / sqrt(I1)
  
  dec <- interim_analysis(ctx$res, z1_obs = z1,
                          theta_hat = theta[THETA_HAT_FIELDS],
                          verbose = FALSE,
                          override_stop = fut == "nonbinding")
  
  n1 <- n_stage1(); nc1 <- 2L * K1
  if ((dec$decision != "continue" && fut == "binding") ||
      (dec$decision == "efficacy_stop" && fut == "nonbinding"))
    return(data.table(decision = dec$decision,
                      reject = dec$decision == "efficacy_stop",
                      z1 = z1, z2 = NA_real_, cstat = NA_real_,
                      I1 = I1, I2 = NA_real_,
                      k2 = 0L, m2 = 0L, N = n1, n_clus = nc1,
                      cost = trial_cost(n1, nc1, 2L)))
  
  k2 <- as.integer(dec$stage2_design$k2)
  m2 <- as.integer(dec$stage2_design$m2)
  t2 <- as.integer(dec$stage2_design$t2 %||% T2)
  
  cl_new <- if (k2 > 0L)
    unlist(lapply(0:1, function(s) setdiff(which(eff$sqv == s), cl1)[seq_len(k2)]))
  else integer(0)
  
  per_cont <- 2L:(1L + t2)
  cpm2 <- sim_periods_cp(eff, cl1, periods = per_cont,
                         m_by_period = rep(m2, length(per_cont)),
                         delta = delta_true, icc = icc_true, p0 = p0_true)
  if (k2 > 0L)
    cpm2 <- rbind(cpm2, sim_periods_cp(eff, cl_new, periods = new_periods(t2),
                                       m_by_period = rep(m2, length(new_periods(t2))),
                                       delta = delta_true, icc = icc_true,
                                       p0 = p0_true))
  
  cpm    <- rbind(cpm1, cpm2)
  Ufulli <- score_contribs(cpm, psi, A_hat, theta$cac, ngrid)
  
  ## U_2|1 = U_full - U_1, aligned on cluster id. New clusters contribute 0 to U_1.
  U1_pad <- setNames(numeric(length(Ufulli)), names(Ufulli))
  U1_pad[names(U1i)] <- U1i
  U2i <- Ufulli - U1_pad
  
  I2 <- sum(U2i^2)
  z2 <- sum(U2i) / sqrt(I2)
  cstat <- W1 * z1 + W2 * z2
  
  n_tot  <- n1 + n_stage2(k2, m2, t2)
  nc_tot <- 2L * (K1 + k2)
  data.table(decision = dec$decision,
             reject = abs(cstat) > ctx$c2,
             z1 = z1, z2 = z2, cstat = cstat, I1 = I1, I2 = I2,
             k2 = k2, m2 = m2, N = n_tot, n_clus = nc_tot,
             cost = trial_cost(n_tot, nc_tot, n_periods_fun(t2)))
}

run_scenario_exact <- function(delta_true, icc_true, cac_true = CAC_P,
                         p0_true = BASELINE_P, ctx = CTX$nonbinding,
                         nsim = NSIM, fut = "nonbinding") {
  rbindlist(lapply_pb(seq_len(nsim), function(i)
    run_once_exact(delta_true, icc_true, cac_true, p0_true, ctx, fut)))
}

library(pbapply)
library(parallel)

set.seed(1234, kind = "L'Ecuyer-CMRG")   # once, at the top of the script
pboptions(type = "timer")                # bar with ETA (off by default in non-interactive sessions)

run_scenario_exact <- function(delta_true, icc_true, cac_true = CAC_P,
                               p0_true = BASELINE_P, ctx = CTX$nonbinding,
                               nsim = NSIM, fut = "nonbinding",
                               cl = cl) {
  # One independent RNG stream per replicate, drawn from the session seed
  stopifnot(RNGkind()[1] == "L'Ecuyer-CMRG")
  s <- get(".Random.seed", envir = globalenv())
  seeds <- vector("list", nsim)
  for (i in seq_len(nsim)) seeds[[i]] <- s <- nextRNGStream(s)
  
  res <- pblapply(seq_len(nsim), function(i) {
    assign(".Random.seed", seeds[[i]], envir = globalenv())
    run_once_exact(delta_true, icc_true, cac_true, p0_true, ctx, fut)
  }, cl = cl)
  
  # Move the session seed past these streams so the next scenario gets fresh ones
  assign(".Random.seed", nextRNGStream(s), envir = globalenv())
  
  # rbindlist silently drops NULLs from crashed workers, so check first
  bad <- vapply(res, function(r) is.null(r) || inherits(r, "try-error"), logical(1))
  if (any(bad)) stop(sum(bad), " of ", nsim, " replicates failed. First error: ", res[bad][[1]])
  rbindlist(res)
}

# -----------------------------------------------------------------------------
# 7. Run the grid
# -----------------------------------------------------------------------------
## 2 hypotheses x 3 ICC x 3 baselines x 2 futility contexts x NSIM.
## Set NSIM <- 200 first to time it before committing to the full run.

## Swap every cached model builder in CTX for a fresh, empty one
refresh_builders <- function(x) {
  if (!is.list(x) || is.data.frame(x)) return(x)
  if (is.function(x[["model_builder"]])) x[["model_builder"]] <- make_crossover_model_fn()
  for (i in seq_along(x))
    if (is.list(x[[i]]) && !is.data.frame(x[[i]])) x[[i]] <- refresh_builders(x[[i]])
  x
}

NCORES <- 7
## Defined at top level and given `fut` rather than `ctx`, so each worker uses its own CTX
run_once_worker <- function(seed, delta_true, icc_true, cac_true, p0_true, fut) {
  assign(".Random.seed", seed, envir = globalenv())
  run_once_exact(delta_true, icc_true, cac_true, p0_true, CTX[[fut]], fut)
}

Sys.setenv(OMP_NUM_THREADS = "1")   # otherwise each fresh worker process uses every core
CL <- makeCluster(NCORES)
invisible(clusterEvalQ(CL, { library(data.table); library(glmmrBase); library(acrt); NULL }))
clusterExport(CL, setdiff(ls(globalenv()),                                   # functions and constants,
                          c("CL", "design_single", "results_single", "results_binding")))  # minus copies already in CTX
invisible(clusterEvalQ(CL, { CTX <- refresh_builders(CTX); NULL }))

## Check: one replicate on a worker must match the master exactly
args <- list(seed = nextRNGStream(.Random.seed), delta_true = DELTA, icc_true = ICC_P,
             cac_true = CAC_P, p0_true = BASELINE_P, fut = "nonbinding")
old <- .Random.seed
on_master <- do.call(run_once_worker, args)
assign(".Random.seed", old, envir = globalenv())   # leave the master's streams untouched
on_worker <- clusterCall(CL[1], do.call, run_once_worker, args)[[1]]
stopifnot(isTRUE(all.equal(on_master, on_worker)))

run_scenario_exact <- function(delta_true, icc_true, cac_true = CAC_P,
                               p0_true = BASELINE_P, ctx = CTX$nonbinding,  # unused: workers look up CTX[[fut]]
                               nsim = NSIM, fut = "nonbinding", cl = CL) {
  stopifnot(RNGkind()[1] == "L'Ecuyer-CMRG")
  s <- get(".Random.seed", envir = globalenv())
  seeds <- vector("list", nsim)
  for (i in seq_len(nsim)) seeds[[i]] <- s <- nextRNGStream(s)
  
  res <- pblapply(seeds, run_once_worker, delta_true = delta_true, icc_true = icc_true,
                  cac_true = cac_true, p0_true = p0_true, fut = fut, cl = cl)
  
  assign(".Random.seed", nextRNGStream(s), envir = globalenv())
  bad <- vapply(res, function(r) is.null(r) || inherits(r, "try-error"), logical(1))
  if (any(bad)) stop(sum(bad), " of ", nsim, " replicates failed. First error: ", res[bad][[1]])
  rbindlist(res)
}

## ... run the simulations, then stopCluster(CL)

sims <- rbindlist(lapply("nonbinding", function(fut) {
  ctx <- CTX[[fut]]
  SCENARIOS[1:18, {
    cat(sprintf("  [%s] %s, ICC = %.2f, CAC = %.2f, p0 = %.2f\n",
                fut, hypothesis, icc_true, cac_true, p0_true))
    run_scenario_exact(delta_true, icc_true, cac_true, p0_true, ctx, NSIM, fut)
  }, by = .(hypothesis, icc_true, cac_true, p0_true)][, futility := fut][]
}))

saveRDS(sims, "impala_adaptive_results.rds")
stopCluster(cl)

## ---- 
# read back in the results
sims1 <- readRDS("C:/Dropbox/impala_adaptive_results.rds")
sims2 <- readRDS("C:/Dropbox/impala_adaptive_results_ne.rds")


# -----------------------------------------------------------------------------
# 8. Summary
# -----------------------------------------------------------------------------
summarise_sims <- function(x, ref_icc = ICC_P, ref_cac = CAC_P,
                           ref_p0 = BASELINE_P) {
  s <- x[, .(
    p_futility  = mean(decision == "futility_stop"),
    p_efficacy  = mean(decision == "efficacy_stop"),
    p_continue  = mean(decision == "continue"),
    reject_rate = mean(reject),
    mcse        = sqrt(mean(reject) * (1 - mean(reject)) / .N),
    E_N   = mean(N), sd_N = sd(N),
    N_q05 = quantile(N, 0.05), N_q95 = quantile(N, 0.95), N_max = max(N),
    E_cost = mean(cost)
  ), by = .(hypothesis, icc_true, cac_true, p0_true, futility)]
  ref <- s[icc_true == ref_icc & cac_true == ref_cac & p0_true == ref_p0,
           .(ref_N = E_N[1]), by = .(hypothesis, futility)]
  s <- merge(s, ref, by = c("hypothesis", "futility"), all.x = TRUE)
  s[, dEN_pct := 100 * (E_N - ref_N) / ref_N][, ref_N := NULL]
  setorder(s, futility, hypothesis, p0_true, icc_true, cac_true)
  s[]
}

summ1     <- summarise_sims(sims1)
summ2 <- summarise_sims(sims2)


summ1$dgp <- "exponential_decay"
summ2$dgp <- "nested_exchangeable"
cols <- colnames(summ1)[c(1,3,4,5,6,7,8,9,11,18)]
dt <- rbind(summ1[,..cols], summ2[,..cols])

library(flextable)

tab <- as.data.frame(dt)

fmt_ci <- function(p, n = 10000) {
  ci <- binom.test(round(p * n), n, conf.level = 0.95)$conf.int
  sprintf("%.3f (%.3f–%.3f)", p, ci[1], ci[2])
}

tab[5:8] <- lapply(tab[5:8], function(x) {
  vapply(x, fmt_ci, character(1))
})

ft <- flextable(tab)
ft <- fontsize(ft, size = 9, part = "all")
ft <- set_table_properties(ft, layout = "autofit", width = 1)

save_as_docx(
  ft, path = "C:/Dropbox/impala_power_results_2.docx",
  pr_section = officer::prop_section(
    page_size = officer::page_size(orient = "landscape")
  )
)



