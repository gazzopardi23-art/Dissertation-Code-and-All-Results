suppressPackageStartupMessages({
  library(data.table)
  library(lubridate)
  library(extRemes)
})

#Helper functions

# Pretty model names
pretty_model_name <- function(nm) {
  switch(nm,
         stat      = "Model 1: Stationary GEV",
         mu        = "Model 2: \u03bc(t) (location varies)",
         sig       = "Model 3: \u03c3(t) (scale varies)",
         xi        = "Model 4: \u03be(t) (shape varies)",
         mu_sig    = "Model 5: \u03bc(t), \u03c3(t)",
         mu_xi     = "Model 6: \u03bc(t), \u03be(t)",
         sig_xi    = "Model 7: \u03c3(t), \u03be(t)",
         mu_sig_xi = "Model 8: \u03bc(t), \u03c3(t), \u03be(t)",
         nm)
}

# GEV density
dgev_custom <- function(x, loc, scale, shape) {
  if (scale <= 0) stop("Scale must be positive.")
  
  if (abs(shape) < 1e-8) {
    z    <- (x - loc) / scale
    dens <- (1 / scale) * exp(-(z + exp(-z)))
  } else {
    tval <- 1 + shape * (x - loc) / scale
    dens <- rep(0, length(x))
    ok   <- tval > 0
    dens[ok] <- (1 / scale) * tval[ok]^(-1 / shape - 1) *
      exp(-tval[ok]^(-1 / shape))
  }
  dens
}

# Extract fitted parameters at a chosen scaled year value
get_params_at_year_sc <- function(fit, year_sc_val) {
  pars <- fit$results$par
  nms  <- names(pars)
  
  get_first <- function(candidates, default = NA_real_) {
    for (nm in candidates) {
      if (nm %in% nms) return(unname(as.numeric(pars[nm])))
    }
    default
  }
  
  loc0 <- get_first(c("location", "mu0"), default = NA_real_)
  loc1 <- get_first(c("location:year_sc", "mu1"), default = 0)
  
  sc0  <- get_first(c("scale", "sigma0"), default = NA_real_)
  sc1  <- get_first(c("scale:year_sc", "sigma1"), default = 0)
  
  sh0  <- get_first(c("shape", "xi0"), default = NA_real_)
  sh1  <- get_first(c("shape:year_sc", "xi1"), default = 0)
  
  if (is.na(loc0)) stop(paste("Missing location parameter. Available:", paste(nms, collapse = ", ")))
  if (is.na(sc0))  stop(paste("Missing scale parameter. Available:",    paste(nms, collapse = ", ")))
  if (is.na(sh0))  stop(paste("Missing shape parameter. Available:",    paste(nms, collapse = ", ")))
  
  loc   <- loc0 + loc1 * year_sc_val
  scale <- sc0  + sc1  * year_sc_val
  shape <- sh0  + sh1  * year_sc_val
  
  if (scale <= 0) {
    stop(sprintf("Non-positive scale at year_sc=%.4f: %.6f", year_sc_val, scale))
  }
  
  list(location = loc, scale = scale, shape = shape)
}

# Delta method CI for return level at a specific year_sc
delta_rl_ci <- function(fit, rp, year_sc_val, alpha = 0.05) {
  
  theta <- fit$results$par
  p     <- length(theta)
  
  rl_fun <- function(pars) {
    nms <- names(theta)
    names(pars) <- nms
    
    loc0 <- if ("location" %in% nms) pars["location"] else pars["mu0"]
    loc1 <- if ("location:year_sc" %in% nms) pars["location:year_sc"] else
      if ("mu1" %in% nms) pars["mu1"] else 0
    
    sc0  <- if ("scale" %in% nms) pars["scale"] else pars["sigma0"]
    sc1  <- if ("scale:year_sc" %in% nms) pars["scale:year_sc"] else
      if ("sigma1" %in% nms) pars["sigma1"] else 0
    
    sh0  <- if ("shape" %in% nms) pars["shape"] else pars["xi0"]
    sh1  <- if ("shape:year_sc" %in% nms) pars["shape:year_sc"] else
      if ("xi1" %in% nms) pars["xi1"] else 0
    
    mu_t  <- as.numeric(loc0 + loc1 * year_sc_val)
    sig_t <- as.numeric(sc0  + sc1  * year_sc_val)
    xi_t  <- as.numeric(sh0  + sh1  * year_sc_val)
    
    yp <- -log(1 - 1 / rp)
    
    if (abs(xi_t) < 1e-8) {
      zT <- mu_t - sig_t * log(yp)
    } else {
      zT <- mu_t + (sig_t / xi_t) * (yp^(-xi_t) - 1)
    }
    zT
  }
  
  zT_hat <- rl_fun(theta)
  
  V <- tryCatch(parcov.fevd(fit), error = function(e) NULL)
  if (is.null(V) && !is.null(fit$results$hessian)) {
    V <- tryCatch(solve(fit$results$hessian), error = function(e) NULL)
  }
  if (is.null(V)) V <- fit$results$cov
  if (is.null(V) || !is.matrix(V) || any(!is.finite(V))) {
    warning(paste("Could not obtain covariance matrix; returning point estimate only"))
    return(c(estimate = unname(zT_hat), lower = NA_real_, upper = NA_real_, se = NA_real_))
  }
  
  eps  <- 1e-5
  grad <- numeric(p)
  for (i in seq_len(p)) {
    theta_up      <- theta
    theta_dn      <- theta
    theta_up[i]   <- theta[i] + eps
    theta_dn[i]   <- theta[i] - eps
    grad[i]       <- (rl_fun(theta_up) - rl_fun(theta_dn)) / (2 * eps)
  }
  
  var_zT <- as.numeric(t(grad) %*% V %*% grad)
  se_zT  <- sqrt(max(var_zT, 0))
  z_crit <- qnorm(1 - alpha / 2)
  
  c(estimate = unname(zT_hat),
    lower    = unname(zT_hat - z_crit * se_zT),
    upper    = unname(zT_hat + z_crit * se_zT),
    se       = unname(se_zT))
}


# Density panel

plot_density_common_scale <- function(fit, bm, model_key, bandwidth = NULL) {
  
  dens_emp <- if (is.null(bandwidth)) {
    density(bm$Tmax, na.rm = TRUE)
  } else {
    density(bm$Tmax, na.rm = TRUE, bw = bandwidth)
  }
  
  x_grid <- seq(min(dens_emp$x), max(dens_emp$x), length.out = 600)
  
  if (model_key == "stat") {
    params   <- get_params_at_year_sc(fit, year_sc_val = 0)
    dens_mod <- dgev_custom(x_grid,
                            loc   = params$location,
                            scale = params$scale,
                            shape = params$shape)
  } else {
    dens_matrix <- sapply(bm$year_sc, function(ysc) {
      params <- tryCatch(
        get_params_at_year_sc(fit, year_sc_val = ysc),
        error = function(e) NULL
      )
      if (is.null(params)) return(rep(NA_real_, length(x_grid)))
      dgev_custom(x_grid,
                  loc   = params$location,
                  scale = params$scale,
                  shape = params$shape)
    })
    dens_mod <- rowMeans(dens_matrix, na.rm = TRUE)
  }
  
  plot(dens_emp,
       main = "",
       xlab = "Annual Max Temperature (\u00b0C)",
       ylab = "Density",
       lwd  = 1)
  
  lines(x_grid, dens_mod, lty = 2, lwd = 1.6, col = "blue")
  
  legend("topright",
         legend = c("Empirical", "Modeled"),
         lty    = c(1, 2),
         lwd    = c(1, 1.6),
         col    = c("black", "blue"),
         bty    = "n")
}

# Save PP, QQ, and custom density in one figure
save_diag_pp_qq_density <- function(fit, bm, file, model_name, model_key) {
  
  png(file, width = 1200, height = 450, res = 140)
  
  op <- par(no.readonly = TRUE)
  on.exit({ par(op); dev.off() }, add = TRUE)
  
  par(mfrow    = c(1, 3),
      mar      = c(5.8, 4.2, 1.8, 1.0),
      oma      = c(0, 0, 2.0, 0),
      mgp      = c(2.4, 0.8, 0),
      cex.axis = 0.95,
      cex.lab  = 1.05)
  
  plot(fit, type = "probprob", main = "")
  mtext("a) Probability plot", side = 1, line = 4.5, font = 2, adj = 0.5)
  
  plot(fit, type = "qq", main = "")
  mtext("b) Q-Q plot", side = 1, line = 4.5, font = 2, adj = 0.5)
  
  plot_density_common_scale(fit = fit, bm = bm, model_key = model_key)
  mtext("c) Empirical distribution", side = 1, line = 4.5, font = 2, adj = 0.5)
  
  mtext(model_name, side = 3, outer = TRUE, line = 0.4, cex = 1.1, font = 2)
}

# Time-varying return levels plot
save_rl_selected_timevarying <- function(fit, bm, RPs, file, model_name) {
  
  n <- nrow(bm)
  
  rl_list <- lapply(RPs, function(rp) {
    v <- as.numeric(return.level(fit, return.period = rp, CI = FALSE))
    if (length(v) == 1) v <- rep(v, n)
    v
  })
  
  rl_mat         <- do.call(cbind, rl_list)
  colnames(rl_mat) <- paste0(RPs, "-year")
  cols           <- c("black", "red3", "blue3", "darkgreen", "purple")
  
  png(file, width = 1100, height = 650, res = 140)
  op <- par(no.readonly = TRUE)
  on.exit({ par(op); dev.off() }, add = TRUE)
  
  par(mar      = c(4.5, 4.5, 3.8, 8.0),
      mgp      = c(2.6, 0.8, 0),
      cex.axis = 0.95,
      cex.lab  = 1.10)
  
  ylim_rng <- range(c(bm$Tmax, rl_mat), na.rm = TRUE)
  
  plot(bm$year, bm$Tmax, type = "l", col = "grey70",
       xlab = "Year", ylab = "Annual Max Temperature (\u00b0C)",
       main = paste0(model_name, " \u2014 Time-varying return levels (2\u201350 years)"),
       ylim = ylim_rng)
  
  for (j in seq_along(RPs)) lines(bm$year, rl_mat[, j], col = cols[j], lwd = 2)
  
  legend("topright", legend = paste0(RPs, "-year"),
         col = cols, lty = 1, lwd = 2, bty = "n",
         inset = c(-0.25, 0), xpd = TRUE)
}

# Stationary return level curve
save_rl_curve_stationary_only <- function(fit, file, model_name) {
  png(file, width = 900, height = 650, res = 140)
  op <- par(no.readonly = TRUE)
  on.exit({ par(op); dev.off() }, add = TRUE)
  
  par(mar      = c(4.5, 4.5, 3.8, 1.2),
      mgp      = c(2.6, 0.8, 0),
      cex.axis = 0.95,
      cex.lab  = 1.10)
  
  plot(fit, type = "rl", main = paste0(model_name, " \u2014 Return level curve"))
}


# 1) Data


dt <- as.data.table(Daily_Temperature_Data_De_Bilt_Station_1950_2024_to_use)
dt[, date  := ymd(sprintf("%08d", as.integer(`YYYYMMDD`)))]
dt[, tx_c  := (`TX - max temperature`) / 10]
dt[, year  := year(date)]


# 2) Block maxima + scaled time


bm <- dt[year %between% c(1950, 2024) & is.finite(tx_c),
         .(Tmax = max(tx_c, na.rm = TRUE)),
         by = year][order(year)]

yrs_ref <- 1950:2024
mu_y    <- mean(yrs_ref)
sd_y    <- sd(yrs_ref)
bm[, year_sc := (year - mu_y) / sd_y]


# 3) GEV fits


model_specs <- list(
  stat      = list(location = ~ 1,        scale = ~ 1,        shape = ~ 1),
  mu        = list(location = ~ year_sc,  scale = ~ 1,        shape = ~ 1),
  sig       = list(location = ~ 1,        scale = ~ year_sc,  shape = ~ 1),
  xi        = list(location = ~ 1,        scale = ~ 1,        shape = ~ year_sc),
  mu_sig    = list(location = ~ year_sc,  scale = ~ year_sc,  shape = ~ 1),
  mu_xi     = list(location = ~ year_sc,  scale = ~ 1,        shape = ~ year_sc),
  sig_xi    = list(location = ~ 1,        scale = ~ year_sc,  shape = ~ year_sc),
  mu_sig_xi = list(location = ~ year_sc,  scale = ~ year_sc,  shape = ~ year_sc)
)

fit_list <- lapply(names(model_specs), function(nm) {
  sp <- model_specs[[nm]]
  fevd(Tmax, data = bm, type = "GEV", method = "MLE",
       location = sp$location,
       scale    = sp$scale,
       shape    = sp$shape)
})
names(fit_list) <- names(model_specs)

fit_stat   <- fit_list[["stat"]]
fit_mu     <- fit_list[["mu"]]
fit_mu_sig <- fit_list[["mu_sig"]]


# 4) AIC/BIC and LRTs


aic_tab <- do.call(rbind, lapply(names(fit_list), function(nm) {
  f <- fit_list[[nm]]
  data.frame(Model = nm,
             k     = length(f$results$par),
             NLL   = f$results$value)
}))
aic_tab$AIC <- 2 * aic_tab$k + 2 * aic_tab$NLL
aic_tab$BIC <- log(nrow(bm)) * aic_tab$k + 2 * aic_tab$NLL
aic_tab     <- aic_tab[order(aic_tab$AIC), ]

lrt_pair <- function(nm0, nm1, label) {
  f0 <- fit_list[[nm0]]; f1 <- fit_list[[nm1]]
  df <- length(f1$results$par) - length(f0$results$par)
  D  <- 2 * (f0$results$value - f1$results$value)
  p  <- 1 - pchisq(D, df = df)
  data.frame(Comparison = label, D = D, df = df, p_value = p)
}

lrt_tab <- rbind(
  lrt_pair("stat",   "mu",        "stat vs mu"),
  lrt_pair("stat",   "sig",       "stat vs sig"),
  lrt_pair("stat",   "xi",        "stat vs xi"),
  lrt_pair("mu",     "mu_sig",    "mu vs mu_sig"),
  lrt_pair("mu",     "mu_xi",     "mu vs mu_xi"),
  lrt_pair("sig",    "mu_sig",    "sig vs mu_sig"),
  lrt_pair("sig",    "sig_xi",    "sig vs sig_xi"),
  lrt_pair("mu_sig", "mu_sig_xi", "mu_sig vs mu_sig_xi"),
  lrt_pair("mu_xi",  "mu_sig_xi", "mu_xi vs mu_sig_xi"),
  lrt_pair("sig_xi", "mu_sig_xi", "sig_xi vs mu_sig_xi")
)


# 5) Parameter estimates + CIs


grab_ci_param <- function(fit) {
  est <- fit$results$par
  out <- try(extRemes::ci(fit, type = "parameter"), silent = TRUE)
  
  if (!inherits(out, "try-error")) {
    M      <- as.matrix(out)
    cn     <- tolower(colnames(M))
    low_i  <- which(grepl("low|2.5",  cn))[1]
    up_i   <- which(grepl("upp|97.5", cn))[1]
    data.frame(Parameter = names(est),
               Estimate  = as.numeric(est),
               lower     = as.numeric(M[, low_i]),
               upper     = as.numeric(M[, up_i]))
  } else {
    covm <- fit$results$cov
    if (!is.null(covm) && length(diag(covm)) == length(est)) {
      se <- sqrt(diag(covm)); z <- qnorm(0.975)
      data.frame(Parameter = names(est),
                 Estimate  = as.numeric(est),
                 lower     = as.numeric(est - z * se),
                 upper     = as.numeric(est + z * se))
    } else {
      data.frame(Parameter = names(est),
                 Estimate  = as.numeric(est),
                 lower     = NA_real_,
                 upper     = NA_real_)
    }
  }
}

par_list <- lapply(fit_list, grab_ci_param)


# 6) Return levels at 2024 with delta method CI for z_50


RPs <- c(2, 5, 10, 20, 50)

year_2024_sc <- (2024 - mu_y) / sd_y

rl_wide_tab <- do.call(rbind, lapply(names(fit_list), function(nm) {
  f <- fit_list[[nm]]
  
  vals <- sapply(RPs, function(rp) {
    delta_rl_ci(f, rp, year_2024_sc)["estimate"]
  })
  
  z50_ci <- delta_rl_ci(f, rp = 50, year_sc_val = year_2024_sc)
  
  data.frame(Model     = nm,
             z2        = vals[1],
             z5        = vals[2],
             z10       = vals[3],
             z20       = vals[4],
             z50       = z50_ci["estimate"],
             z50_lower = z50_ci["lower"],
             z50_upper = z50_ci["upper"],
             z50_se    = z50_ci["se"],
             row.names = NULL)
}))

rl_wide_tab <- rl_wide_tab[order(rl_wide_tab$Model), ]


# 7) Save outputs


base_out <- "bm_output_debilt_tmax_FIXED"
dir.create(base_out, showWarnings = FALSE)

fwrite(bm[, .(year, Tmax)], file.path(base_out, "DeBilt_BM_Tmax_1950_2024.csv"))
fwrite(aic_tab,             file.path(base_out, "GEV_AIC_BIC_all_models.csv"))
fwrite(lrt_tab,             file.path(base_out, "LRT_key_comparisons.csv"))

for (nm in names(par_list)) {
  fwrite(par_list[[nm]], file.path(base_out, sprintf("Params_%s.csv", nm)))
}

fwrite(rl_wide_tab, file.path(base_out, "ReturnLevels_WIDE_all_models_2024.csv"))


# 8) Save figures


png(file.path(base_out, "BM_series.png"), width = 900, height = 600, res = 140)
plot(bm$year, bm$Tmax, pch = 16,
     xlab = "Year", ylab = "Annual Max Temperature (\u00b0C)",
     main = "De Bilt Annual Maximum Temperature (1950\u20132024)")
abline(lm(Tmax ~ year, bm), lty = 2)
dev.off()

diag_dir <- file.path(base_out, "diagnostics_common_scale")
dir.create(diag_dir, showWarnings = FALSE)

for (nm in names(fit_list)) {
  save_diag_pp_qq_density(
    fit        = fit_list[[nm]],
    bm         = bm,
    file       = file.path(diag_dir, sprintf("Diagnostics_COMMONSCALE_%s.png", nm)),
    model_name = pretty_model_name(nm),
    model_key  = nm
  )
}

rl_time_dir <- file.path(base_out, "returnlevels_selected_2_50_timevarying")
dir.create(rl_time_dir, showWarnings = FALSE)

for (nm in names(fit_list)) {
  save_rl_selected_timevarying(
    fit        = fit_list[[nm]],
    bm         = bm,
    RPs        = RPs,
    file       = file.path(rl_time_dir, sprintf("RL_SELECTED_2_50_%s.png", nm)),
    model_name = pretty_model_name(nm)
  )
}

rl_curve_dir <- file.path(base_out, "returnlevel_curve_stationary")
dir.create(rl_curve_dir, showWarnings = FALSE)

save_rl_curve_stationary_only(
  fit        = fit_stat,
  file       = file.path(rl_curve_dir, "RLcurve_stat.png"),
  model_name = pretty_model_name("stat")
)


# 9) Console summary


cat("\n=== AIC/BIC (sorted by AIC) ===\n")
print(aic_tab, row.names = FALSE)

cat("\n=== LRT key comparisons ===\n")
print(lrt_tab, row.names = FALSE)

cat("\n=== Return levels at 2024 with delta method 95% CI for z_50 ===\n")
print(rl_wide_tab, row.names = FALSE)