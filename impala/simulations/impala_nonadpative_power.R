## require Packages

library(lme4)
library(dplyr)
library(tidyr)
library(glmmrBase)
library(glmmTMB)
library(statmod)

set.seed(564354)

# this code calculates the power for IMPALA 2 using multiple methods. The code below will reproduce the calculations
# from www.clustertrial.app. Note that the "design effect" method replicates the approach in the Shiny app.

### FUNCTIONS ####

expit <- function(x) 1 / (1 + exp(-x))

# turn 1 trial to four non-overlapping
mult_data <- function(df){
  dfo <- df
  for(i in 2:4){
    df2 <- dfo
    df2$t <- df2$t + max(df$t)
    df2$cl <- df2$cl + max(df$cl)
    df <- rbind(df,df2)
  }
  return(df)
}

# convert parameter values - marginal
to_logit_scale <- function(baseline, eff_size){
  b0 <- log(baseline/(1-baseline))
  m1 <- baseline*eff_size
  b1 <- log(m1/(1-m1)) - b0
  return(c(b0,b1))
}

## GET THE GLMM PARAMETERS FROM THE SPECIFIED CORRELATION PARAMETERS
## conditional logit parameters from specified marginal probabilities

solve_parameters_xsect_binomial <- function(p0, p1, icc, max_iter = 100, tol = 1e-6,
                                            n_quad = 30) {
  
  # --- Gauss-Hermite quadrature for N(0,1) integration ---
  gh <- statmod::gauss.quad(n_quad, kind = "hermite")
  nodes   <- gh$nodes * sqrt(2)
  weights <- gh$weights / sqrt(pi)
  
  compute_moments <- function(eta, sigma_c) {
    p   <- plogis(eta + sigma_c * nodes)
    mu  <- sum(weights * p)
    mu2 <- sum(weights * p^2)
    list(mean = mu, var_p = mu2 - mu^2)
  }
  
  compute_icc <- function(var_p, mu) {
    var_p / (mu * (1 - mu))
  }
  
  # --- Initial guesses (Zeger-Liang-style attenuation) ---
  pi2_3 <- pi^2 / 3
  cc    <- 0.588
  
  sigma_c_init <- sqrt(icc * pi2_3 / (1 - icc + 0.01))
  sigma_c_init <- max(0.1, min(2.0, sigma_c_init))
  
  atten      <- sqrt(1 + cc^2 * sigma_c_init^2)
  beta0_init <- log(p0 / (1 - p0)) * atten
  beta1_init <- log(p1 / (1 - p1)) * atten - beta0_init
  
  params <- c(beta0_init, beta1_init, log(sigma_c_init))
  
  beta_max        <- 10
  log_sigma_range <- c(-4, 4)
  eps    <- 1e-5
  lambda <- 0.1
  
  # --- Levenberg-Marquardt loop ---
  for (iter in seq_len(max_iter)) {
    sigma_c <- exp(params[3])
    
    m0 <- compute_moments(params[1], sigma_c)
    m1 <- compute_moments(params[1] + params[2], sigma_c)
    
    curr_p0  <- m0$mean
    curr_p1  <- m1$mean
    curr_icc <- compute_icc(m0$var_p, m0$mean)
    
    resid  <- c(p0 - curr_p0, p1 - curr_p1, icc - curr_icc)
    sum_sq <- sum(resid^2)
    
    if (sqrt(sum_sq / 3) < tol) {
      return(list(beta0 = params[1], beta1 = params[2],
                  sigma_c = sigma_c, sigma_t = 0,
                  converged = TRUE, iterations = iter))
    }
    
    # Numerical Jacobian (3 x 3)
    J <- matrix(0, 3, 3)
    for (j in 1:3) {
      pp <- params
      h  <- max(eps, abs(params[j]) * eps)
      pp[j] <- pp[j] + h
      
      sc_p <- exp(pp[3])
      m0_p <- compute_moments(pp[1], sc_p)
      m1_p <- compute_moments(pp[1] + pp[2], sc_p)
      icc_p <- compute_icc(m0_p$var_p, m0_p$mean)
      
      J[1, j] <- (m0_p$mean - curr_p0) / h
      J[2, j] <- (m1_p$mean - curr_p1) / h
      J[3, j] <- (icc_p - curr_icc) / h
    }
    
    # LM update
    JTJ <- crossprod(J)
    JTr <- crossprod(J, resid)
    diag(JTJ) <- diag(JTJ) + lambda
    delta <- solve(JTJ, JTr)
    
    # Trust region clamp
    sn <- sqrt(sum(delta^2))
    if (sn > 2) delta <- delta * 2 / sn
    
    # Backtracking line search
    alpha    <- 1
    improved <- FALSE
    
    for (ls in 1:15) {
      pn <- params + alpha * as.numeric(delta)
      pn[1] <- max(-beta_max, min(beta_max, pn[1]))
      pn[2] <- max(-beta_max, min(beta_max, pn[2]))
      pn[3] <- max(log_sigma_range[1], min(log_sigma_range[2], pn[3]))
      
      sc_n   <- exp(pn[3])
      m0_n   <- compute_moments(pn[1], sc_n)
      m1_n   <- compute_moments(pn[1] + pn[2], sc_n)
      icc_n  <- compute_icc(m0_n$var_p, m0_n$mean)
      new_sq <- (p0 - m0_n$mean)^2 + (p1 - m1_n$mean)^2 + (icc - icc_n)^2
      
      if (new_sq < sum_sq) {
        params   <- pn
        improved <- TRUE
        lambda   <- max(1e-8, lambda * 0.7)
        break
      }
      alpha <- alpha * 0.5
    }
    
    if (!improved) lambda <- min(1e6, lambda * 3)
  }
  
  sigma_c <- exp(params[3])
  list(beta0 = params[1], beta1 = params[2],
       sigma_c = sigma_c, sigma_t = 0,
       converged = FALSE, iterations = max_iter)
}

# for the design effect we need the correlation matrix
build_gee_covariance_matrix <- function(data, icc, iac, cac, mean_n, mu0, mu1,
                                        correlation_structure = c("nested_exchangeable", "exponential_decay")) {
  
  correlation_structure <- match.arg(correlation_structure)
  
  n_obs <- nrow(data)
  cluster <- data[, 1]
  period  <- data[, 2]
  
  mu_bar   <- (mu0 + mu1) / 2
  totalvar <- mu_bar * (1 - mu_bar)
  
  sig2CP <- icc * totalvar
  sig2E  <- (1 - iac) * (totalvar - sig2CP)
  sig2   <- sig2E / mean_n
  sigindiv <- if (iac > 0 && iac < 1) {
    sig2E * iac / ((1 - iac) * mean_n)
  } else {
    0
  }
  
  # Edge cases for IAC
  if (iac == 0) {
    sig2E    <- totalvar - sig2CP
    sig2     <- sig2E / mean_n
    sigindiv <- 0
  } else if (iac >= 1 - 1e-10) {
    sig2E    <- 0
    sig2     <- 0
    sigindiv <- (totalvar - sig2CP) / mean_n
  }
  
  r <- cac
  
  # Build covariance matrix for cluster-period means
  V <- matrix(0, n_obs, n_obs)
  
  for (i in seq_len(n_obs)) {
    for (j in i:n_obs) {
      if (cluster[i] != cluster[j]) {
        cov_ij <- 0
      } else {
        lag <- abs(period[i] - period[j])
        if (correlation_structure == "nested_exchangeable") {
          cov_ij <- if (lag == 0) {
            sigindiv + sig2 + sig2CP
          } else {
            sigindiv + r * sig2CP
          }
        } else {
          cov_ij <- if (lag == 0) {
            sigindiv + sig2 + sig2CP
          } else {
            sigindiv + sig2CP * r^lag
          }
        }
      }
      V[i, j] <- cov_ij
      V[j, i] <- cov_ij
    }
  }
  
  V
}

# calculate the design effect standard error for the power
gee_variance <- function(X, Sigma_GEE, cv = 0, idx = NULL) {
  
  n_obs <- nrow(X)
  
  # CV correction to diagonal
  if (cv > 0) {
    diag(Sigma_GEE) <- diag(Sigma_GEE) * (1 + cv^2)
  }
  
  # Check positive definiteness
  chol_Sigma <- tryCatch(chol(Sigma_GEE), error = function(e) {
    stop("GEE covariance matrix not positive definite")
  })
  
  # GLS information matrix: X' Sigma^{-1} X
  # Using the Cholesky factor to solve
  L_inv_X <- forwardsolve(t(chol_Sigma), X)
  M_gee <- crossprod(L_inv_X)
  Minv <- solve(M_gee)
  
  # Extract variance for target coefficient
  if (is.null(idx)) idx <- ncol(X)
  bvar <- Minv[idx, idx]
  
  if (is.nan(bvar) || bvar <= 0) {
    stop("Invalid GEE variance estimate")
  }
  
  list(
    se   = sqrt(bvar),
    Minv = Minv,
    M    = M_gee
  )
}


##############################
# Calculate all the powers          #
##############################



power <- function(p0, p1, K, 
                  n_periods, n, icc, ar, 
                  cov_dgp = c("exponential_decay", "nested_exchangeable"),
                  niter = 1000, cl = NULL){
  # generate dummy trial data
  
  # n_periods is number of trial periods (crossovers), either 2 or 4
  # K is number of clusters per arm
  # n is cluster-period size
  cat("\nGenerate data")
  gen_data <- function(n_periods, K, n){
    df <- nelder(as.formula(paste0("~ cl(",K*2,") * t(",n_periods,")")))
    df$int <- 0
    df[df$cl <= K & df$t %in% c(1,3), 'int'] <- 1
    df[df$cl > K & df$t %in% c(2,4), 'int'] <- 1
    df$n <- n
    df
  }
  
  ind_data <- function(data){
    df[rep(1:nrow(data),data$n),]
  }
  
  df <- gen_data(n_periods,K,n)
  cat("\nSolve for covariance parameters")
  pars <-  solve_parameters_xsect_binomial(p0,p1,icc)
  
  cat("\nForm model")
  if(cov_dgp == "exponential_decay"){
    cl_mod <- Model$new(
      ~ int + factor(t) + (1|gr(cl)*ar0(t)),
      data = df, 
      family = binomial(),
      mean = c(pars$beta0,pars$beta1,rep(0,n_periods - 1)),
      covariance = c(pars$sigma_c,ar),
      trials = rep(n,nrow(df))
    )
  } else if(cov_dgp == "nested_exchangeable"){
    
    cl_mod <- Model$new(
      ~ int + factor(t) + (1|gr(cl)) + (1|gr(cl,t)),
      data = df, 
      family = binomial(),
      mean = c(pars$beta0,pars$beta1,rep(0,n_periods - 1)),
      covariance = c(pars$sigma_c*ar,pars$sigma_c*(1-ar)),
      trials = rep(n,nrow(df))
    )
  }
  
  M <- solve(cl_mod$information_matrix(average = FALSE))
  
  cat("\nCalculate power")
  V_gee <- build_gee_covariance_matrix(df,icc,0,ar,n,p0,p1,cov_dgp)
  se_gee <- gee_variance(cl_mod$mean$X,V_gee,cv=0,idx = 2)
  pow_de <- pnorm(abs((p1 - p0)/se_gee$se) - qnorm(1- 0.05/2))
  
  # GLS
  pow_gls <- pnorm(abs(pars$beta1) / sqrt(M[2,2]) - qnorm(0.975))
  
  # t-test
  pow_gls_t <- pt(abs(pars$beta1) / sqrt(M[2,2]) - qt(0.975, df = K*2-1-n_periods), df = K*2-1-n_periods)
  
  # kenward-roger
  Mkr <- cl_mod$small_sample_correction(type = "KR")
  pow_kr <- pt(abs(pars$beta1) / sqrt(Mkr$vcov_beta[2,2]) - qt(0.975, df = Mkr$dof[2]), df = Mkr$dof[2])
  
  # simulation-based
  cat("\nRun simulation...\n")
  res <- data.frame(b = rep(NA, niter), se = rep(NA, niter))
  # outcomes 
  outcomes <- matrix(NA, nrow = nrow(df), ncol = niter)
  for(i in 1:niter) outcomes[,i] <- cl_mod$sim_data()
  df$t_f <- numFactor(df$t)
  fit_ar <- function(i, df, outcomes){
    df$y <- outcomes[,i]
    fit <- glmmTMB(cbind(y, n - y) ~ int + factor(t) + ar1(t_f + 0 | cl), family= binomial(), data = df)
    sfit <- summary(fit)
    return(c(sfit$coefficients$cond[2,1], sfit$coefficients$cond[2,2]))
  }
  
  if(!require(glmmTMB))parallel::clusterEvalQ(cl, .libPaths("C:/R/lib")) # needed to set path on my laptop
  parallel::clusterEvalQ(cl, require(glmmTMB))
  parallel::clusterExport(cl, c("outcomes", "df"), envir = environment())
  result <- pbapply::pbsapply(1:niter, function(x) fit_ar(x,df,outcomes), cl = cl)
  res$b <- result[1,]
  res$se <- result[2,]
  res$t <- res$b/res$se
  res$p <- 2*(1 - pt(abs(res$t), df = K*2-1-n_periods))
  pow_sim <- mean(res$p < 0.05)
  cover_sim <- mean(res$b + qt(0.975, df = K*2-1-n_periods)*res$se > pars$beta1 &
                      res$b - qt(0.975, df = K*2-1-n_periods)*res$se < pars$beta1)
  return(list(pow_de = pow_de, pow_gls = pow_gls, pow_gls_t = pow_gls_t, 
              pow_kr = pow_kr, pow_sim = pow_sim,
              cover_sim = cover_sim, res = res))
}


scenarios <- expand.grid(baseline = c(0.12, 0.15, 0.18),
                         effect = c(-0.02, -0.03, -0.04),
                         icc = c(0.05, 0.10, 0.15),
                         ar = c(0.8, 0.9),
                         n = c(75, 150),
                         dgp = c("exponential_decay"))

scenarios <- rbind(scenarios,expand.grid(baseline = c(0.15),
                                         effect = c(-0.03),
                                         icc = c(0.05, 0.10, 0.15),
                                         ar = c(0.8, 0.9),
                                         n = c(75),
                                         dgp = c("nested_exchangeable")))



cl <- parallel::makeCluster(7)
range1 <- 1:38
range2 <- 39:114
results <- list()
for(i in range2){
  cat("\nIter ",i," of ",nrow(scenarios))
  p1 <- power(scenarios$baseline[i],
                scenarios$baseline[i] + scenarios$effect[i],
                15, 4, scenarios$n[i],
                scenarios$icc[i],
                scenarios$ar[i],
                as.character(scenarios$dgp[i]),niter = 10000,cl)
  results[[i - min(range2) + 1]] <- list(params = scenarios[i,], power = p1)
  saveRDS(results,"impala_power_results_part2.rds")
}

parallel::stopCluster(cl)


## combine results

df <- readRDS("C:/Dropbox/impala_power_results_part1.rds")
df2 <- readRDS("C:/Dropbox/impala_power_results_part2.rds")

results <- append(df, df2)

# install.packages(c("flextable", "officer"))
library(flextable)

fmt_ci <- function(p, n = 10000) {
  ci <- binom.test(round(p * n), n, conf.level = 0.95)$conf.int
  sprintf("%.3f (%.3f–%.3f)", p, ci[1], ci[2])
}

tab <- do.call(rbind, lapply(results, function(x) {
  p <- x$power
  data.frame(
    x$params,
    DE  = sprintf("%.3f", p$pow_de),
    GLS = sprintf("%.3f", p$pow_gls),
    `GLS (t)` = sprintf("%.3f", p$pow_gls_t),
    KR  = sprintf("%.3f", p$pow_kr),
    `Sim. power (95% CI)` = fmt_ci(p$pow_sim),
    `Coverage (95% CI)`   = fmt_ci(p$cover_sim),
    check.names = FALSE
  )
}))

ft <- flextable(tab)
ft <- fontsize(ft, size = 9, part = "all")
ft <- set_table_properties(ft, layout = "autofit", width = 1)

save_as_docx(
  ft, path = "C:/Dropbox/impala_power_results.docx",
  pr_section = officer::prop_section(
    page_size = officer::page_size(orient = "landscape")
  )
)
