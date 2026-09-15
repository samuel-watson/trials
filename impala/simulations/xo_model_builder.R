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

convert_to_glmm_scale <- function(family = c("gaussian", "binomial", "poisson"),
                                 baseline = 0,
                                 delta,
                                 icc,
                                 cac = NULL) {
  
  family <- match.arg(family)
  
  if (family == "gaussian") {
    
    beta <- c(baseline, delta)
    
    cov_pars <- if (!is.null(cac)) {
      c(icc * cac, icc * (1 - cac))
    } else {
      icc
    }
    
    return(list(
      beta = beta,
      cov_pars = cov_pars,
      var_par = 1 - icc,
      family = gaussian(),
      b1_original = delta,
      b1_link = delta,
      link = "identity"
    ))
    
  } else if (family == "binomial") {
    
    p0 <- baseline
    p1 <- baseline + delta
    
    if (p0 <= 0 || p0 >= 1)
      stop("baseline must lie in (0,1).")
    
    if (p1 <= 0 || p1 >= 1)
      stop("baseline + delta must lie in (0,1).")
    
    ## Existing approximate marginal-ICC -> latent variance mapping
    A <- (icc / (1 - icc)) / (p0 * (1 - p0))
    
    ## Exact marginal mean calibration conditional on A
    lk <- calibrate_link(p0, delta, A)
    
    beta <- c(lk$eta0, lk$b1)
    
    cov_pars <- if (!is.null(cac)) {
      c(A * cac, A * (1 - cac))
    } else {
      A
    }
    
    return(list(
      beta = beta,
      cov_pars = cov_pars,
      tau_sq = A,
      var_par = 1,
      family = binomial(),
      b1_original = delta,
      b1_link = lk$b1,
      link = "logit"
    ))
    
  } else if (family == "poisson") {
    
    b0 <- log(baseline)
    rate1 <- baseline + delta
    
    if (rate1 <= 0)
      stop("baseline + delta must be positive.")
    
    b1 <- log(rate1) - b0
    tau_sq <- (icc / (1 - icc)) / baseline
    
    cov_pars <- if (!is.null(cac)) {
      c(tau_sq * cac, tau_sq * (1 - cac))
    } else {
      tau_sq
    }
    
    return(list(
      beta = c(b0, b1),
      cov_pars = cov_pars,
      tau_sq = tau_sq,
      var_par = 1,
      family = poisson(),
      b1_original = delta,
      b1_link = b1,
      link = "log"
    ))
  }
}

crossover_model_builder <- function(design_params, fixed_params) {
  cache <- fixed_params$cache
  k1 <- fixed_params$k1
  m11 <- fixed_params$m11
  m12 <- fixed_params$m12
  k2 <- design_params$k2
  m2 <- design_params$m2
  m2b <- design_params$m2b %||% design_params$m2
  icc <- fixed_params$icc
  cac <- fixed_params$cac
  delta <- fixed_params$delta
  cache <- fixed_params$cache
  t2 <- design_params$t2
  # Family and baseline
  family <- fixed_params$family %||% "gaussian"
  baseline <- fixed_params$baseline %||% 0
  
  # Convert parameters to appropriate scale
  glm_params <- convert_to_glmm_scale(
    family = family,
    baseline = baseline,
    delta = delta,
    icc = icc,
    cac = cac
  )
  # Stage 1: crossover (2 periods)
  df1 <- glmmrBase::nelder(as.formula(paste0("~ (int(2) > cl(", k1, ")) * t(2)")))
  df1$int <- df1$int - 1
  df1$t <- df1$t - 1
  # Treatment: arm 1 gets trt in period 2, arm 2 gets trt in period 1
  df1$trt <- ifelse((df1$int == 0 & df1$t == 1) | (df1$int == 1 & df1$t == 0), 1, 0)
  df1$n <- m11
  df1[df1$t == 1, "n"] <- m12
  if(t2 == 1){
    # Stage 2 period (period 3)
    df2_base <- glmmrBase::nelder(as.formula(paste0("~ int(2) > cl(", k1, ")")))
    df2_base$int <- df2_base$int - 1
    df2_base$t <- 2
    df2_base$trt <- df2_base$int  # Same as their original arm
    df2_base$n <- m2
  } else {
    df2_base <- glmmrBase::nelder(as.formula(paste0("~ (int(2) > cl(", k1, ")) * t(2)")))
    df2_base$int <- df2_base$int - 1
    df2_base$trt <- ifelse((df2_base$int == 0 & df2_base$t == 2) | (df2_base$int == 1 & df2_base$t == 1), 1, 0)
    df2_base$n <- m2
    df2_base$t <- df2_base$t + 1
  }
  
  if (k2 > 0) {
    # New clusters get both periods
    df2_new <- glmmrBase::nelder(as.formula(paste0("~ (int(2) > cl(", k2, ")) * t(2)")))
    df2_new$int <- df2_new$int - 1
    df2_new$trt <- ifelse((df2_new$int == 0 & df2_new$t == 2) | (df2_new$int == 1 & df2_new$t == 1), 1, 0)
    df2_new$n <- m2
    df2_new$t <- df2_new$t + 1
    df2_new$cl <- df2_new$cl + max(df1$cl)
    df <- rbind(df1, df2_base, df2_new)
  } else {
    df <- rbind(df1, df2_base)
  }
  # Number of time periods
  n_periods <- length(unique(df$t))
  
  # Mean parameters: intercept + (n_periods - 1) time effects + treatment + extra 0 for random
  # glm_params$beta = c(intercept, treatment)
  # Need: c(intercept, time2, time3, ..., treatment, 0)
  beta_full <- c(glm_params$beta[1],           # Intercept
                 rep(0, n_periods - 1),         # Time effects
                 glm_params$beta[2]         )   # Treatment
  
  # Cache key includes family
  cache_key <- paste(k1, k2, t2, family, sep = "_")
  if (is.null(cache$mod) || cache$cache_key != cache_key) {
    
    if (family == "gaussian") {
      cache$mod <- glmmrBase::Model$new(
        ~ factor(t) + trt + (1|gr(cl)*ar0(t)),
        data = df,
        family = glm_params$family,
        mean = beta_full,
        covariance = glm_params$cov_pars,
        weights = df$n,
        var_par = glm_params$var_par
      )
    } else if (family == "binomial") {
      print(c(glm_params$cov_pars[1]/cac, cac))
      print(beta_full)
      cache$mod <- glmmrBase::Model$new(
        ~ factor(t) + trt + (1|gr(cl)*ar0(t)),
        data = df,
        family = glm_params$family,
        mean = beta_full,
        covariance =  c(glm_params$cov_pars[1]/cac, cac),
        trials = df$n
      )
    } else if (family == "poisson") {
      cache$mod <- glmmrBase::Model$new(
        ~ factor(t) + trt + (1|gr(cl)*ar0(t)),
        data = df,
        family = glm_params$family,
        mean = beta_full,
        covariance = glm_params$cov_pars,
        offset = log(df$n)
      )
    }
    
    cache$cache_key <- cache_key
    
  } else {
    # Update existing model
    if (family == "gaussian") {
      glmmrBase:::Model__set_weights(cache$mod$.__enclos_env__$private$ptr, df$n)
    } else if (family == "binomial") {
      glmmrBase:::Model__set_trials(cache$mod$.__enclos_env__$private$ptr, as.integer(df$n))
    } else if (family == "poisson") {
      glmmrBase:::Model__set_offset(cache$mod$.__enclos_env__$private$ptr, log(df$n))
    }
    cache$mod$update_parameters(cov.pars = c(glm_params$cov_pars[1]/cac, cac), mean.pars = beta_full)
  }
  n1 <- nrow(df1)
  n2 <- nrow(df) - n1
  X <- cache$mod$mean$X
  D <- cache$mod$covariance$D
  Z <- cache$mod$covariance$Z
  Z_sp <- Matrix::Matrix(Z, sparse = TRUE)
  D_sp <- Matrix::Matrix(D, sparse = TRUE)
  S <- Matrix::Diagonal(x = 1/cache$mod$w_matrix()) + Z_sp %*% Matrix::tcrossprod(D_sp, Z_sp)
  V <- as.matrix(S)
  
  idx1 <- 1:n1
  idx2 <- (n1 + 1):nrow(df)
  
  # Treatment column - find it by name
  j_trt <- which(colnames(X) == "trt")
  if (length(j_trt) == 0) {
    # Fallback: last column before random effects
    j_trt <- ncol(X)
  }
  eff_decomp <- efficient_score_decomposition(X, V, idx1, idx2, j = j_trt)
  # Degrees of freedom
  # Cluster-periods: stage 1 has 2*k1 clusters * 2 periods, stage 2 adds k1 + 2*k2 cluster-periods
  n_cluster_periods_s1 <- 2 * k1 * 2  # 2 arms * k1 clusters * 2 periods
  n_cluster_periods_full <- n_cluster_periods_s1 + t2 * 2 * (k1 + k2)  # + stage 2
  # Fixed effects: intercept + (n_periods - 1) time effects + treatment
  n_fixed <- 1 + (n_periods - 1) + 1
  
  df_s1 <- n_cluster_periods_s1 - (1 + 1 + 1)  # Only 2 periods in stage 1, so fewer time effects
  
  df_full <- n_cluster_periods_full - n_fixed
  list(
    I1_eff = eff_decomp$I1_eff,
    I2_eff = eff_decomp$I2_eff,
    I_eff = eff_decomp$I_eff,
    w1 = sqrt(eff_decomp$I1_eff / eff_decomp$I_eff),
    w2 = sqrt(eff_decomp$I2_eff / eff_decomp$I_eff),
    b1 = glm_params$b1_link %||% glm_params$beta[2],
    b1_original = glm_params$b1_original,
    n1 = n1,
    n2 = n2,
    m11 = m11,
    m12 = m12,
    m2 = m2,
    m2b = m2b,
    k1 = k1,
    k2 = k2,
    t2 = t2,
    df_s1 = df_s1,
    df_s2 = t2 * 2 * (k1 + k2) - t2 - 1,
    df_full = df_full,
    family = family,
    link = glm_params$link
  )
}

# Also need factory function for caching
make_crossover_model_fn <- function() {
  cache <- new.env()
  cache$mod <- NULL
  cache$cache_key <- NULL
  
  function(design_params, fixed_params) {
    fixed_params$cache <- cache
    crossover_model_builder(design_params, fixed_params)
  }
}

crossover_crt <- function(icc, cac = 0.8, delta,
                          k1, m11, m12,
                          k2 = 0:4, m2 = NULL, t2 = 1:2,
                          rho = 30, rho_t = 10,
                          family = "gaussian", baseline = NULL) {
  
  spec <- crt_design_spec(
    stage1_params = c("k1", "m11", "m12"),
    stage2_params = c("k2", "m2", "t2"),
    resources = list(
      n_s1 = ~ n_arms * k1 * (m11 + m12),        # 2 periods in stage 1
      n_s2 = ~ n_arms * k1 * m2 * t2 + n_arms * k2 * m2 * t2,  # existing + new clusters
      clusters_s1 = ~ n_arms * k1,
      clusters_s2 = ~ n_arms * k2,
      t_s1 = ~ 2,
      t_s2 = ~ t2
    ),
    cost_structure = list(
      weights = c(n = 1, clusters = "rho", t = "rho_t"),
      stage2_resources = c("n_s2", "clusters_s2", "t_s2")
    ),
    model_builder = crossover_model_builder,
    n_arms = 2,
    design_type = "crossover"
  )
  structure(
    list(
      spec = spec,
      fixed_params = list(
        icc = icc,
        cac = cac,
        delta = delta,
        family = family,
        baseline = baseline %||% 0,
        rho_t = rho_t
      ),
      stage1_grid = expand.grid(k1 = k1, m11 = m11, m12 = m12),
      stage2_grid_fn = function(s1) expand.grid(t2 = t2, k2 = k2, m2 = m2),
      rho = rho
    ),
    class = c("crossover_crt", "crt_design")
  )
}

