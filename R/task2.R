# ============================================================================
# MAST32001 - Task 2: Bayesian mixture modelling with MCMC
# Gibbs sampling, simulation validation and Old Faithful analysis.
#
# Log:
# 16.10.2026: folders are not being created in the correct path; add a
#             helper function to locate the project root and set the
#             working directory; finished the gibbs sampler; part B in progress.
# 17.10.2026: complete part B analysis and update helper functions.
#             generated csv and PNG outputs are not included.
#             part C in progress.
# 18.10.2026: use global theme palettes; standardize section titles
#             and descriptive figure/table output names.
#             reorder colors and enlarge text for report figures.
# 19.10.2026: run two manually initialized chains, use PNG null devices,
#             and add numerical summaries for the complete results report.
# 20.10.2026：resize charts, relabel sections, and reuse project_setup.R
#             instead of first three sections; save all four additional summary
#             tables in the workspace; handle missing or non-finite ESS values
#             in diagnostic screening.
# 22.10.2026: append the team's code after section 6.7 and save
#             the added results in task2_sim.RData; remove section 8.4.
# 24.10.2026: resize charts.
# ============================================================================

## ---- 0. shared project setup ----------------------------------------------
root <- normalizePath("...")
source(file.path(root, file.path("R", "project_setup.R")))
ensure_packages(c("posterior"))

## ---- 1. Part A: Gibbs sampler implementation ------------------------------
# Full-conditional Gibbs sampler for the 2-component Gaussian mixture.
rgibbs_mix <- function(y,
                       n_iter = 5000,
                       burn = 1000,
                       m0 = 0,
                       s0 = 10,
                       a = 1,
                       b = 1,
                       alpha0 = 2,
                       beta0 = 1,
                       init = NULL,
                       seed = NULL) {
  if (!is.null(seed))
    set.seed(seed)
  y <- as.numeric(y)
  n <- length(y)
  scalar <- function(x)
    is.numeric(x) &&
    length(x) == 1 &&
    is.finite(x)
  if (n < 3 ||
      any(!is.finite(y)) || sd(y) <= 0)
    stop("Invalid data")
  if (!scalar(n_iter) || !scalar(burn) || n_iter != floor(n_iter) ||
      burn != floor(burn) || burn < 0 || n_iter <= burn) {
    stop("Require integer n_iter > burn >= 0")
  }
  if (!all(vapply(list(m0, s0, a, b, alpha0, beta0), scalar, logical(1))) ||
      any(c(s0, a, b, alpha0, beta0) <= 0))
    stop("Invalid prior parameters")
  started <- proc.time()[["elapsed"]]
  mu1 <- mu2 <- pi <- s2 <- numeric(n_iter)
  z1 <- integer(n_iter)                       # n1 at each iteration
  
  if (is.null(init)) {
    # k-means start (10 restarts)
    km <- kmeans(y, centers = 2, nstart = 10)
    mu1_c <- km$centers[1]
    mu2_c <- km$centers[2]
    pi_c  <- km$size[1] / n
    s2_c  <- mean((y - km$centers[km$cluster])^2)
    if (s2_c <= 0)
      s2_c <- var(y) / 2
  } else {
    mu1_c <- init$mu1
    mu2_c <- init$mu2
    pi_c  <- init$pi
    s2_c  <- init$s2
  }
  
  initial <- c(mu1_c, mu2_c, pi_c, s2_c)
  if (length(initial) != 4 || any(!is.finite(initial)) ||
      pi_c <= 0 ||
      pi_c >= 1 || s2_c <= 0)
    stop("Invalid initial values")
  init_used <- setNames(as.numeric(initial), c("mu1", "mu2", "pi", "s2"))
  
  for (it in seq_len(n_iter)) {
    # 1. Allocations
    log1 <- log(pi_c)   + dnorm(y, mu1_c, sqrt(s2_c), log = TRUE)
    log2 <- log1p(-pi_c) + dnorm(y, mu2_c, sqrt(s2_c), log = TRUE)
    p1   <- plogis(log1 - log2)
    zi   <- rbinom(n, 1, p1)
    n1   <- sum(zi == 1)
    n2 <- n - n1
    
    # 2. Component means
    ybar1 <- if (n1 > 0)
      mean(y[zi == 1])
    else
      0
    ybar2 <- if (n2 > 0)
      mean(y[zi == 0])
    else
      0
    prec1 <- 1 / s0^2 + n1 / s2_c
    mu1_c <- rnorm(1, (m0 / s0^2 + n1 * ybar1 / s2_c) / prec1, sqrt(1 /
                                                                      prec1))
    prec2 <- 1 / s0^2 + n2 / s2_c
    mu2_c <- rnorm(1, (m0 / s0^2 + n2 * ybar2 / s2_c) / prec2, sqrt(1 /
                                                                      prec2))
    
    # 3. Mixing proportion
    pi_c <- rbeta(1, a + n1, b + n2)
    
    # 4. Variance: mean priors are independent N(m0,s0^2).
    # IG density is proportional to v^(-alpha-1)*exp(-beta/v).
    # There is no extra K/2 in the conditional shape.
    resid <- sum(ifelse(zi == 1, (y - mu1_c)^2, (y - mu2_c)^2))
    s2_c  <- 1 / rgamma(1, shape = alpha0 + n / 2, rate = beta0 + resid /
                          2)
    if (any(!is.finite(c(mu1_c, mu2_c, pi_c, s2_c))) ||
        s2_c <= 0 ||
        pi_c <= 0 ||
        pi_c >= 1)
      stop("Invalid sampled state at iteration ", it)
    
    mu1[it] <- mu1_c
    mu2[it] <- mu2_c
    pi[it]  <- pi_c
    s2[it]  <- s2_c
    z1[it]  <- n1
  }
  
  out <- data.frame(
    mu1 = mu1,
    mu2 = mu2,
    pi = pi,
    s2 = s2,
    n1 = z1
  )
  
  list(
    samples = out[seq.int(burn + 1, n_iter), , drop = FALSE],
    all = out,
    init_used = init_used,
    n = n,
    n_iter = n_iter,
    burn = burn,
    seed = seed,
    elapsed_seconds = proc.time()[["elapsed"]] - started,
    prior = list(
      m0 = m0,
      s0 = s0,
      a = a,
      b = b,
      alpha0 = alpha0,
      beta0 = beta0
    )
  )
}

## ---- 2. Part B: simulation settings and parameter recovery ----------------
# ---- 2.1 Separated and overlapping simulation settings ---------------------
n <- 500
set.seed(20261016)
z_true <- rbinom(n, 1, 0.5)
y_sep <- ifelse(z_true == 1, rnorm(n, -2, 1), rnorm(n, 2, 1))    # mu = +/- 2
y_ov  <- ifelse(z_true == 1, rnorm(n, -0.5, 1), rnorm(n, 0.5, 1))  # mu = +/- 0.5

# ---- 2.2 Run chains from different initial values --------------------------
inits <- list(
  # swapped labels
  swap   = list(
    mu1 = 2,
    mu2 = -2,
    pi = 0.9,
    s2 = 1
  ),
  # spread means
  spread = list(
    mu1 = -4,
    mu2 = 4,
    pi = 0.5,
    s2 = 4
  )
)
run_chains <- function(y,
                       inits,
                       n_iter = 5000,
                       burn = 1000,
                       seed_base = 261024100) {
  ans <- lapply(seq_along(inits), function(i) {
    message("Chain ", i, ": ", names(inits)[i])
    rgibbs_mix(
      y,
      n_iter = n_iter,
      burn = burn,
      init = inits[[i]],
      seed = seed_base + i
    )
  })
  setNames(ans, names(inits))
}
message("Separated simulation")
chains_sep <- run_chains(y_sep,
                         inits,
                         n_iter = 20000,
                         burn = 5000,
                         seed_base = 261024100)
message("Overlapping simulation")
chains_ov  <- run_chains(y_ov,
                         inits,
                         n_iter = 60000,
                         burn = 15000,
                         seed_base = 261024100)

# ---- 2.3 Posterior relabeling after sampling -------------------------------
relabel_chain <- function(ch) {
  s <- ch$samples
  flip <- s$mu1 > s$mu2
  data.frame(
    mu_lo = pmin(s$mu1, s$mu2),
    mu_hi = pmax(s$mu1, s$mu2),
    pi_lo = ifelse(flip, 1 - s$pi, s$pi),
    s2 = s$s2,
    n_lo = ifelse(flip, ch$n - s$n1, s$n1)
  )
}
pool_relabelled <- function(chains) {
  do.call(rbind, lapply(chains, function(ch) {
    relabel_chain(ch)[, c("mu_lo", "mu_hi", "pi_lo", "s2")]
  }))
}
rl_sep <- pool_relabelled(chains_sep)
rl_ov  <- pool_relabelled(chains_ov)

# ---- 2.4 Posterior allocation uncertainty ----------------------------------
alloc_post <- function(ch, y, thin = 1) {
  if (!scalar_numeric(thin) ||
      thin < 1 || thin != floor(thin))
    stop("Invalid thinning interval")
  s <- relabel_chain(ch)
  keep <- seq(1, nrow(s), by = thin)
  s <- s[keep, , drop = FALSE]
  vapply(y, function(yi) {
    log1 <- log(s$pi_lo) + dnorm(yi, s$mu_lo, sqrt(s$s2), log = TRUE)
    log2 <- log1p(-s$pi_lo) + dnorm(yi, s$mu_hi, sqrt(s$s2), log = TRUE)
    mean(plogis(log1 - log2))
  }, numeric(1))
}
scalar_numeric <- function(x)
  is.numeric(x) && length(x) == 1 && is.finite(x)
alloc_chains <- function(chains, y) {
  P <- vapply(chains, alloc_post, numeric(length(y)), y = y)
  w <- vapply(chains, function(ch)
    nrow(ch$samples), integer(1))
  as.numeric(P %*% (w / sum(w)))
}
allocation_table <- function(y, p) {
  safe <- pmin(1 - 1e-15, pmax(1e-15, p))
  data.frame(
    observation = seq_along(y),
    y = y,
    p_low = p,
    ambiguous_02_08 = p > 0.2 & p < 0.8,
    entropy = -safe * log(safe) - (1 - safe) * log1p(-safe)
  )
}
p_sep <- alloc_chains(chains_sep, y_sep)
p_ov  <- alloc_chains(chains_ov, y_ov)
alloc_sep <- allocation_table(y_sep, p_sep)
alloc_ov <- allocation_table(y_ov, p_ov)
write.csv(alloc_sep,
          file.path(table_dir, "task2_separated_allocations.csv"),
          row.names = FALSE)
write.csv(alloc_ov,
          file.path(table_dir, "task2_overlap_allocations.csv"),
          row.names = FALSE)
amb_sep <- mean(p_sep > 0.2 & p_sep < 0.8)
amb_ov  <- mean(p_ov  > 0.2 & p_ov  < 0.8)

# ---- 2.5 Simulation parameter recovery summaries ---------------------------
mk_rec <- function(truth, mu_lo, mu_hi, pi_lo, s2) {
  data.frame(
    Parameter = c("mu1 (lower)", "mu2 (upper)", "pi", "sigma2"),
    Truth = truth,
    Median = c(median(mu_lo), median(mu_hi), median(pi_lo), median(s2)),
    lo = c(
      quantile(mu_lo, 0.025),
      quantile(mu_hi, 0.025),
      quantile(pi_lo, 0.025),
      quantile(s2, 0.025)
    ),
    hi = c(
      quantile(mu_lo, 0.975),
      quantile(mu_hi, 0.975),
      quantile(pi_lo, 0.975),
      quantile(s2, 0.975)
    )
  )
}
tab_sep <- mk_rec(c(-2, 2, 0.5, 1),
                  rl_sep$mu_lo,
                  rl_sep$mu_hi,
                  rl_sep$pi_lo,
                  rl_sep$s2)
tab_ov  <- mk_rec(c(-0.5, 0.5, 0.5, 1),
                  rl_ov$mu_lo,
                  rl_ov$mu_hi,
                  rl_ov$pi_lo,
                  rl_ov$s2)
write.csv(
  tab_sep,
  file.path(table_dir, "task2_separated_parameter_summary.csv"),
  row.names = FALSE
)
write.csv(tab_ov,
          file.path(table_dir, "task2_overlap_parameter_summary.csv"),
          row.names = FALSE)

## ---- 3. Part B: trace plots before relabelling ----------------------------
# ---- 3.1 Prepare simulation trace data -------------------------------------
trace_frame <- function(chains, setting) {
  do.call(rbind, lapply(seq_along(chains), function(i) {
    ch <- chains[[i]]
    s <- ch$samples
    # Reduce only displayed points; all draws are used for inference.
    keep <- unique(round(seq(1, nrow(s), length.out = min(1500, nrow(
      s
    )))))
    data.frame(
      setting = setting,
      chain = names(chains)[i],
      iter = ch$burn + keep,
      mu1 = s$mu1[keep],
      mu2 = s$mu2[keep]
    )
  }))
}
tr <- rbind(trace_frame(chains_sep, "separated"),
            trace_frame(chains_ov, "overlap")) |>
  mutate(setting = factor(
    setting,
    levels = c("separated", "overlap"),
    labels = c("Separated components", "Overlapping components")
  ))
tr_long <- tidyr::pivot_longer(tr, c(mu1, mu2), names_to = "param", values_to = "value")

# ---- 3.2 Plot means with their original component labels -------------------
p_trace_sep <- ggplot(
  filter(tr_long, setting == "Separated components"),
  aes(
    iter,
    value,
    colour = param,
    linetype = chain,
    group = interaction(param, chain)
  )
) +
  geom_line(linewidth = 0.3, alpha = 0.85) +
  scale_colour_discrete(
    name = NULL,
    limits = c("mu1", "mu2"),
    labels = c(mu1 = expression(mu[1]), mu2 = expression(mu[2]))
  ) +
  scale_linetype_discrete(name = "Chain", limits = c("swap", "spread")) +
  labs(x = "Iteration (post burn-in)", y = "MCMC draw") +
  theme(legend.position = "none")

p_trace_ov <- ggplot(
  filter(tr_long, setting == "Overlapping components"),
  aes(
    iter,
    value,
    colour = param,
    linetype = chain,
    group = interaction(param, chain)
  )
) +
  geom_line(linewidth = 0.3, alpha = 0.85) +
  scale_colour_discrete(
    name = NULL,
    limits = c("mu1", "mu2"),
    labels = c(mu1 = expression(mu[1]), mu2 = expression(mu[2]))
  ) +
  scale_linetype_discrete(name = "Chain", limits = c("swap", "spread")) +
  labs(x = "Iteration (post burn-in)", y = "MCMC draw") +
  theme(legend.position = "none")

fig_simulation_raw_traces <- ggarrange(
  p_trace_sep,
  p_trace_ov,
  ncol = 2,
  labels = c("a", "b"),
  font.label = list(size = 12, face = "bold"),
  label.x = 0,
  label.y = 1,
  hjust = 0,
  vjust = 1,
  common.legend = TRUE,
  legend = "top",
  align = "v"
)
save_sci(
  fig_simulation_raw_traces,
  "task2_simulation_raw_traces",
  width = "double",
  height_mm = 75
)

## ---- 4. Part B: MCMC convergence and efficiency diagnostics ---------------
# ---- 4.1 Rank-normalised Rhat, ESS and MCSE --------------------------------
diag_details <- function(chains, ordered = TRUE) {
  # Retain chain identities when calculating rank-normalised split Rhat and ESS.
  mats <- lapply(chains, function(ch) {
    if (ordered)
      as.matrix(relabel_chain(ch)[, c("mu_lo", "mu_hi", "pi_lo", "s2")])
    else
      as.matrix(ch$samples[, c("mu1", "mu2", "pi", "s2")])
  })
  sizes <- vapply(mats, nrow, integer(1))
  if (length(unique(sizes)) != 1)
    stop("Diagnostic chain lengths must agree")
  arr <- array(NA_real_,
               c(sizes[1], length(mats), 4),
               dimnames = list(NULL, names(chains), colnames(mats[[1]])))
  for (i in seq_along(mats))
    arr[, i, ] <- mats[[i]]
  d <- as.data.frame(
    posterior::summarise_draws(
      posterior::as_draws_array(arr),
      "mean",
      "median",
      "sd",
      "rhat",
      "ess_bulk",
      "ess_tail",
      "mcse_mean"
    )
  )
  names(d) <- c("Parameter",
                "Mean",
                "Median",
                "SD",
                "Rhat",
                "ESS",
                "ESS_tail",
                "MCSE")
  d
}
diag_summary <- function(chains,
                         ordered = TRUE,
                         details = NULL) {
  d <- if (is.null(details))
    diag_details(chains, ordered)
  else
    details
  # Preserve the original wide table format. ESS_mu1/ESS_mu2 are low/high when ordered.
  c(
    ESS_mu1 = d$ESS[1],
    ESS_mu2 = d$ESS[2],
    ESS_pi = d$ESS[3],
    ESS_s2 = d$ESS[4],
    Rhat_max = max(d$Rhat)
  )
}
# ---- 4.2 Simulation diagnostics with ordered and original labels -----------
diag_sep <- diag_details(chains_sep)
diag_ov <- diag_details(chains_ov)
diag_sep_raw <- diag_details(chains_sep, FALSE)
diag_ov_raw <- diag_details(chains_ov, FALSE)
diag_tab <- rbind(
  Separated = diag_summary(chains_sep, details = diag_sep),
  Overlap = diag_summary(chains_ov, details = diag_ov)
)
diag_raw_tab <- rbind(
  Separated = diag_summary(chains_sep, FALSE, diag_sep_raw),
  Overlap = diag_summary(chains_ov, FALSE, diag_ov_raw)
)
diag_detailed <- rbind(cbind(Setting = "Separated", diag_sep),
                       cbind(Setting = "Overlap", diag_ov))
diag_raw_detailed <- rbind(cbind(Setting = "Separated", diag_sep_raw),
                           cbind(Setting = "Overlap", diag_ov_raw))
write.csv(diag_tab,
          file.path(table_dir, "task2_ordered_diagnostics_overview.csv"))
write.csv(diag_raw_tab,
          file.path(table_dir, "task2_raw_label_diagnostics_overview.csv"))
write.csv(diag_detailed,
          file.path(table_dir, "task2_ordered_diagnostics.csv"),
          row.names = FALSE)
write.csv(
  diag_raw_detailed,
  file.path(table_dir, "task2_raw_label_diagnostics.csv"),
  row.names = FALSE
)

# ---- 4.3 Diagnostic screening thresholds -----------------------------------
diagnostic_screen <- function(d, setting) {
  good <- all(is.finite(d$Rhat), is.finite(d$ESS), is.finite(d$ESS_tail)) &&
    max(d$Rhat) <= 1.01 &&
    min(d$ESS) >= 400 && min(d$ESS_tail) >= 400
  if (!good)
    warning(setting,
            ": inspect traces and increase iterations before reporting final results")
  data.frame(
    Setting = setting,
    Rhat_max = max(d$Rhat),
    ESS_min = min(d$ESS),
    ESS_tail_min = min(d$ESS_tail),
    Screen_passed = good
  )
}
screen_tab <- rbind(diagnostic_screen(diag_sep, "Separated"),
                    diagnostic_screen(diag_ov, "Overlap"))

## ---- 5. Part C: label switching and posterior summaries -------------------
# ---- 5.1 Label changes and allocation ambiguity ----------------------------
switches <- function(ch)
  sum(diff(sign(ch$samples$mu1 - ch$samples$mu2)) != 0)
sw_sep <- sapply(chains_sep, switches)
sw_ov  <- sapply(chains_ov, switches)
cat(sprintf(
  "Label switches (per chain): separated %s, overlap %s\n",
  paste(sw_sep, collapse = ","),
  paste(sw_ov, collapse = ",")
))
cat(
  sprintf(
    "Ambiguous allocations (0.2-0.8): separated %.1f%%, overlap %.1f%%\n",
    100 * amb_sep,
    100 * amb_ov
  )
)

# ---- 5.2 Compare summaries with original and ordered labels ----------------
label_summary <- function(chains, setting) {
  do.call(rbind, lapply(seq_along(chains), function(i) {
    ch <- chains[[i]]
    raw <- ch$samples
    ord <- relabel_chain(ch)
    data.frame(
      Setting = setting,
      Chain = names(chains)[i],
      Switches = switches(ch),
      Fraction_mu1_lower = mean(raw$mu1 < raw$mu2),
      Raw_mu1 = median(raw$mu1),
      Raw_mu2 = median(raw$mu2),
      Raw_pi = median(raw$pi),
      Ordered_mu_lo = median(ord$mu_lo),
      Ordered_mu_hi = median(ord$mu_hi),
      Ordered_pi_lo = median(ord$pi_lo)
    )
  }))
}
label_tab <- rbind(label_summary(chains_sep, "Separated"),
                   label_summary(chains_ov, "Overlap"))
# ---- 5.3 Posterior quantities invariant to component labels ----------------
invariant_summary <- function(chains, setting) {
  s <- do.call(rbind, lapply(chains, function(ch)
    ch$samples))
  x <- data.frame(
    mixture_mean = s$pi * s$mu1 + (1 - s$pi) * s$mu2,
    mixture_variance = s$s2 + s$pi * (1 - s$pi) * (s$mu1 - s$mu2)^2,
    separation = abs(s$mu1 - s$mu2),
    s2 = s$s2
  )
  data.frame(
    Setting = setting,
    Quantity = names(x),
    Median = vapply(x, median, numeric(1)),
    lo = vapply(x, quantile, numeric(1), probs = 0.025),
    hi = vapply(x, quantile, numeric(1), probs = 0.975),
    row.names = NULL
  )
}
invariant_tab <- rbind(
  invariant_summary(chains_sep, "Separated"),
  invariant_summary(chains_ov, "Overlap")
)

# ---- 5.4 Simulated data: generating and fitted densities -------------------
mixture_curve <- function(chains, y) {
  s <- do.call(rbind, lapply(chains, function(ch)
    ch$samples))
  xs <- seq(min(y) - 0.15 * diff(range(y)),
            max(y) + 0.15 * diff(range(y)),
            length.out = 250)
  data.frame(x = xs, fit = vapply(xs, function(x)
    mean(
      s$pi * dnorm(x, s$mu1, sqrt(s$s2)) + (1 - s$pi) * dnorm(x, s$mu2, sqrt(s$s2))
    ), numeric(1)))
}
density_sep <- mixture_curve(chains_sep, y_sep)
density_ov <- mixture_curve(chains_ov, y_ov)
simulation_density_plot <- function(y, d, means) {
  d$truth <- 0.5 * dnorm(d$x, means[1], 1) + 0.5 * dnorm(d$x, means[2], 1)
  ggplot(data.frame(y = y), aes(y)) +
    geom_histogram(
      aes(y = after_stat(density), fill = "Observed"),
      bins = 30,
      alpha = 0.25,
      colour = "white",
      linewidth = 0.2
    ) +
    geom_line(
      data = d,
      aes(x, fit, colour = "Fitted"),
      inherit.aes = FALSE,
      linewidth = 0.6
    ) +
    geom_line(
      data = d,
      aes(x, truth, colour = "Generating"),
      inherit.aes = FALSE,
      linewidth = 0.5,
      linetype = 2
    ) +
    scale_colour_discrete(name = NULL, limits = c("Fitted", "Generating")) +
    scale_fill_discrete(guide = "none") +
    labs(x = "Simulated observation", y = "Density")
}
fig_simulation_density <- ggarrange(
  simulation_density_plot(y_sep, density_sep, c(-2, 2)),
  simulation_density_plot(y_ov, density_ov, c(-0.5, 0.5)),
  ncol = 2,
  labels = c("a", "b"),
  font.label = list(size = 12, face = "bold"),
  label.x = 0,
  label.y = 1,
  hjust = 0,
  vjust = 1,
  common.legend = TRUE,
  legend = "top"
)
save_sci(
  fig_simulation_density,
  "task2_simulation_density",
  width = "double",
  height_mm = 75
)

# ---- 5.5 Simulation traces after posterior relabelling ---------------------
ordered_trace_frame <- function(chains, setting) {
  do.call(rbind, lapply(seq_along(chains), function(i) {
    ch <- chains[[i]]
    s <- relabel_chain(ch)
    keep <- unique(round(seq(1, nrow(s), length.out = min(1500, nrow(
      s
    )))))
    do.call(rbind, lapply(c("mu_lo", "mu_hi", "pi_lo", "s2"), function(p)
      data.frame(
        setting = setting,
        chain = names(chains)[i],
        param = p,
        iter = ch$burn + keep,
        value = s[[p]][keep]
      )))
  }))
}
tr_ordered <- rbind(
  ordered_trace_frame(chains_sep, "Separated"),
  ordered_trace_frame(chains_ov, "Overlap")
)
tr_ordered$setting <- factor(tr_ordered$setting, levels = c("Separated", "Overlap"))
fig_simulation_ordered_traces <- ggplot(tr_ordered, aes(iter, value, colour = chain)) +
  geom_line(linewidth = 0.25, alpha = 0.85) + facet_grid(param ~ setting, scales = "free") +
  scale_colour_discrete(name = "Chain", limits = c("swap", "spread")) +
  labs(x = "Iteration (post burn-in)", y = "Ordered parameter")
save_sci(
  fig_simulation_ordered_traces,
  "task2_simulation_ordered_traces",
  width = "double",
  height_mm = 120
)

## ---- 6. Part D: Old Faithful analysis -------------------------------------
# ---- 6.1 Faithful data and chain initialisation ----------------------------
data(faithful, package = "datasets")
y_faith <- faithful$eruptions
q_f <- as.numeric(quantile(y_faith, c(0.2, 0.8)))
inits_f <- list(
  swap = list(
    mu1 = q_f[2],
    mu2 = q_f[1],
    pi = 0.8,
    s2 = var(y_faith)
  ),
  spread = list(
    mu1 = min(y_faith) - sd(y_faith),
    mu2 = max(y_faith) + sd(y_faith),
    pi = 0.5,
    s2 = 2 * var(y_faith)
  )
)
# ---- 6.2 Sample the Faithful posterior and calculate diagnostics -----------
message("Old Faithful")
chains_f <- run_chains(
  y_faith,
  inits_f,
  n_iter = 20000,
  burn = 5000,
  seed_base = 261024300
)
chain_f <- chains_f[[1]]   # First raw chain retained for compatibility.
rl_f <- pool_relabelled(chains_f)
diag_f <- diag_details(chains_f)
diag_f_raw <- diag_details(chains_f, FALSE)
ess_f <- min(diag_f$ESS)
p_faith_alloc <- alloc_chains(chains_f, y_faith)
amb_faith <- mean(p_faith_alloc > 0.2 & p_faith_alloc < 0.8)
alloc_f <- allocation_table(y_faith, p_faith_alloc)
write.csv(alloc_f,
          file.path(table_dir, "task2_faithful_allocations.csv"),
          row.names = FALSE)

# ---- 6.3 Combine diagnostics and label summaries across datasets -----------
diag_tab <- rbind(diag_tab, Faithful = diag_summary(chains_f, details = diag_f))
diag_raw_tab <- rbind(diag_raw_tab, Faithful = diag_summary(chains_f, FALSE, diag_f_raw))
diag_detailed <- rbind(diag_detailed, cbind(Setting = "Faithful", diag_f))
diag_raw_detailed <- rbind(diag_raw_detailed, cbind(Setting = "Faithful", diag_f_raw))
screen_tab <- rbind(screen_tab, diagnostic_screen(diag_f, "Faithful"))
label_tab <- rbind(label_tab, label_summary(chains_f, "Faithful"))
invariant_tab <- rbind(invariant_tab, invariant_summary(chains_f, "Faithful"))
write.csv(diag_tab,
          file.path(table_dir, "task2_ordered_diagnostics_overview.csv"))
write.csv(diag_raw_tab,
          file.path(table_dir, "task2_raw_label_diagnostics_overview.csv"))
write.csv(diag_detailed,
          file.path(table_dir, "task2_ordered_diagnostics.csv"),
          row.names = FALSE)
write.csv(
  diag_raw_detailed,
  file.path(table_dir, "task2_raw_label_diagnostics.csv"),
  row.names = FALSE
)
write.csv(screen_tab,
          file.path(table_dir, "task2_diagnostic_screen.csv"),
          row.names = FALSE)
write.csv(label_tab,
          file.path(table_dir, "task2_label_switching.csv"),
          row.names = FALSE)
write.csv(
  invariant_tab,
  file.path(table_dir, "task2_label_invariant_summary.csv"),
  row.names = FALSE
)

# ---- 6.4 Faithful posterior parameter summary ------------------------------
tab_f <- data.frame(
  Parameter = c("mu (short)", "mu (long)", "pi (short)", "sigma2"),
  Median = c(
    median(rl_f$mu_lo),
    median(rl_f$mu_hi),
    median(rl_f$pi_lo),
    median(rl_f$s2)
  ),
  lo = c(
    quantile(rl_f$mu_lo, 0.025),
    quantile(rl_f$mu_hi, 0.025),
    quantile(rl_f$pi_lo, 0.025),
    quantile(rl_f$s2, 0.025)
  ),
  hi = c(
    quantile(rl_f$mu_lo, 0.975),
    quantile(rl_f$mu_hi, 0.975),
    quantile(rl_f$pi_lo, 0.975),
    quantile(rl_f$s2, 0.975)
  ),
  check.names = FALSE
)
write.csv(tab_f,
          file.path(table_dir, "task2_faithful_parameter_summary.csv"),
          row.names = FALSE)

# ---- 6.5 Observed Faithful data and fitted mixture density -----------------
s <- do.call(rbind, lapply(chains_f, function(ch)
  ch$samples))
xs <- seq(min(y_faith), max(y_faith), length.out = 400)
f_pred <- sapply(xs, function(x)
  mean(s$pi * dnorm(x, s$mu1, sqrt(s$s2)) + (1 - s$pi) * dnorm(x, s$mu2, sqrt(s$s2))))
fig_faithful_density <- ggplot(data.frame(y = y_faith), aes(y)) +
  geom_histogram(
    aes(y = after_stat(density), fill = "Observed"),
    bins = 30,
    alpha = 0.25,
    colour = "white",
    linewidth = 0.2
  ) +
  geom_line(data = data.frame(x = xs, y = f_pred),
            aes(x, y, colour = "Fitted"),
            linewidth = 0.6) +
  scale_colour_discrete(limits = "Fitted", guide = "none") +
  scale_fill_discrete(guide = "none") +
  labs(x = "Eruption duration (min)", y = "Density")
save_sci(
  fig_faithful_density,
  "task2_faithful_density",
  width = "double",
  height_mm = 75
)

# ---- 6.6 Posterior allocation uncertainty across datasets ------------------
allocation_plot <- function(a, xlab) {
  ggplot(a, aes(y, p_low)) +
    geom_hline(
      yintercept = c(0.2, 0.8),
      linetype = 3,
      colour = "grey60",
      linewidth = 0.25
    ) +
    geom_point(aes(colour = "Allocation"),
               size = 1.1,
               alpha = 0.8) +
    scale_colour_discrete(guide = "none") +
    scale_y_continuous(limits = c(0, 1)) + labs(x = xlab, y = "Pr(lower component | data)")
}
fig_allocation_probabilities <- ggarrange(
  allocation_plot(alloc_sep, "Separated simulation"),
  allocation_plot(alloc_ov, "Overlapping simulation"),
  allocation_plot(alloc_f, "Eruption duration (min)"),
  ncol = 3,
  labels = c("a", "b", "c"),
  font.label = list(size = 12, face = "bold"),
  label.x = 0,
  label.y = 1,
  hjust = 0,
  vjust = 1,
  align = "v"
)
save_sci(
  fig_allocation_probabilities,
  "task2_allocation_probabilities",
  width = "double",
  height_mm = 75
)

# ---- 6.7 Faithful traces after posterior relabeling ------------------------
tr_f <- ordered_trace_frame(chains_f, "Faithful")
fig_faithful_ordered_traces <- ggplot(tr_f, aes(iter, value, colour = chain)) +
  geom_line(linewidth = 0.25, alpha = 0.85) + facet_wrap(~ param, scales = "free_y", ncol = 2) +
  scale_colour_discrete(name = "Chain", limits = c("swap", "spread")) +
  labs(x = "Iteration (post burn-in)", y = "Ordered parameter")
save_sci(
  fig_faithful_ordered_traces,
  "task2_faithful_ordered_traces",
  width = "double",
  height_mm = 120
)
invisible(gc())

## ---- 7. Additions: Part C Relabeling Critique & Part D Faithful Follow-up --

# ---- 7.1 Part C: 2D Joint Posterior (mu1 vs mu2) & Empty-Component Effect ---
joint_trace_df <- function(chains, setting) {
  do.call(rbind, lapply(seq_along(chains), function(i) {
    ch <- chains[[i]]
    s  <- ch$samples
    keep <- unique(round(seq(1, nrow(s), length.out = min(2000, nrow(
      s
    )))))
    data.frame(
      setting = setting,
      chain   = names(chains)[i],
      mu1     = s$mu1[keep],
      mu2     = s$mu2[keep],
      empty   = (s$n1[keep] == 0 | s$n1[keep] == ch$n)
    )
  }))
}

df_joint <- rbind(
  joint_trace_df(chains_sep, "Separated components"),
  joint_trace_df(chains_ov, "Overlapping components")
)
df_joint$setting <- factor(df_joint$setting,
                           levels = c("Separated components", "Overlapping components"))

fig_joint_posterior <- ggplot(df_joint, aes(x = mu1, y = mu2)) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = 2,
    colour = "grey40",
    linewidth = 0.4
  ) +
  geom_point(aes(colour = chain, shape = empty),
             alpha = 0.55,
             size = 1.1) +
  facet_wrap( ~ setting, scales = "free") +
  scale_colour_discrete(name = "Chain", limits = c("swap", "spread")) +
  scale_shape_manual(
    name   = "Empty component",
    values = c(`FALSE` = 16, `TRUE` = 4),
    labels = c(`FALSE` = "No (n1, n2 > 0)", `TRUE` = "Yes (prior draw)")
  ) +
  labs(x = expression(mu[1]), y = expression(mu[2]))

save_sci(
  fig_joint_posterior,
  "task2_joint_posterior_mu1_mu2",
  width = "double",
  height_mm = 75
)

# Compare all ordered draws vs. occupancy-conditioned ordered draws (min(n1, n2) >= 5)
relabel_filtered_summary <- function(chains, min_occ = 5) {
  ord_all <- do.call(rbind, lapply(chains, relabel_chain))
  n_obs   <- chains[[1]]$n
  keep    <- ord_all$n_lo >= min_occ &
    (n_obs - ord_all$n_lo) >= min_occ
  ord_sub <- ord_all[keep, , drop = FALSE]
  
  summarize_set <- function(df, label) {
    params <- c("mu_lo", "mu_hi", "pi_lo", "s2")
    data.frame(
      Scheme        = label,
      Retained_pct  = 100 * nrow(df) / nrow(ord_all),
      Parameter     = params,
      Median        = vapply(df[params], median, numeric(1)),
      lo            = vapply(df[params], quantile, numeric(1), probs = 0.025),
      hi            = vapply(df[params], quantile, numeric(1), probs = 0.975),
      row.names     = NULL
    )
  }
  rbind(
    summarize_set(ord_all, "All ordered draws (mu_lo < mu_hi)"),
    summarize_set(
      ord_sub,
      paste0("Occupancy-filtered (min n_k >= ", min_occ, ")")
    )
  )
}

tab_ov_relabel_critique <- relabel_filtered_summary(chains_ov, min_occ = 5)
write.csv(
  tab_ov_relabel_critique,
  file.path(table_dir, "task2_overlap_relabeling_critique.csv"),
  row.names = FALSE
)

# ---- 7.2 Part D: Formal Posterior Predictive Check (PPC) for Old Faithful ---
run_faithful_ppc <- function(rl_draws,
                             y_obs,
                             split_val = 3.0,
                             n_rep = 2000,
                             seed = 20261022) {
  set.seed(seed)
  n <- length(y_obs)
  idx <- sample(seq_len(nrow(rl_draws)),
                size = min(n_rep, nrow(rl_draws)),
                replace = FALSE)
  
  calc_stats <- function(y_vec) {
    s_short <- sd(y_vec[y_vec < split_val])
    s_long  <- sd(y_vec[y_vec >= split_val])
    c(
      sd_short = s_short,
      sd_long = s_long,
      sd_ratio = s_long / s_short
    )
  }
  
  obs_stats <- calc_stats(y_obs)
  rep_stats <- vapply(idx, function(m) {
    z_rep <- rbinom(n, 1, rl_draws$pi_lo[m])
    y_rep <- ifelse(z_rep == 1,
                    rnorm(n, rl_draws$mu_lo[m], sqrt(rl_draws$s2[m])),
                    rnorm(n, rl_draws$mu_hi[m], sqrt(rl_draws$s2[m])))
    calc_stats(y_rep)
  }, numeric(3))
  
  rep_df <- as.data.frame(t(rep_stats))
  summary_df <- data.frame(
    Statistic = c(
      "SD short (< 3 min)",
      "SD long (>= 3 min)",
      "Ratio (SD long / SD short)"
    ),
    Observed  = as.numeric(obs_stats),
    Rep_Median = vapply(rep_df, median, numeric(1)),
    Rep_lo     = vapply(rep_df, quantile, numeric(1), probs = 0.025),
    Rep_hi     = vapply(rep_df, quantile, numeric(1), probs = 0.975),
    PPP_value  = c(
      mean(rep_df$sd_short >= obs_stats["sd_short"]),
      mean(rep_df$sd_long  >= obs_stats["sd_long"]),
      mean(rep_df$sd_ratio >= obs_stats["sd_ratio"])
    ),
    row.names = NULL
  )
  list(summary = summary_df,
       rep_draws = rep_df,
       obs = obs_stats)
}

ppc_f <- run_faithful_ppc(rl_f, y_faith)
write.csv(ppc_f$summary,
          file.path(table_dir, "task2_faithful_ppc_summary.csv"),
          row.names = FALSE)

# ---- 7.3 Part D: Unequal-Variance Gibbs Sampler on Old Faithful -------------
rgibbs_mix_unequal <- function(y,
                               n_iter = 20000,
                               burn = 5000,
                               m0 = 0,
                               s0 = 10,
                               a = 1,
                               b = 1,
                               alpha0 = 2,
                               beta0 = 1,
                               init = NULL,
                               seed = NULL) {
  if (!is.null(seed))
    set.seed(seed)
  y <- as.numeric(y)
  n <- length(y)
  
  mu1 <- mu2 <- pi <- s2_1 <- s2_2 <- numeric(n_iter)
  z1  <- integer(n_iter)
  
  mu1_c  <- init$mu1
  mu2_c  <- init$mu2
  pi_c   <- init$pi
  s2_1_c <- init$s2
  s2_2_c <- init$s2
  
  for (it in seq_len(n_iter)) {
    # 1. Allocations with component-specific variances
    log1 <- log(pi_c)    + dnorm(y, mu1_c, sqrt(s2_1_c), log = TRUE)
    log2 <- log1p(-pi_c) + dnorm(y, mu2_c, sqrt(s2_2_c), log = TRUE)
    zi   <- rbinom(n, 1, plogis(log1 - log2))
    n1   <- sum(zi == 1)
    n2   <- n - n1
    
    # 2. Component means
    ybar1 <- if (n1 > 0)
      mean(y[zi == 1])
    else
      0
    ybar2 <- if (n2 > 0)
      mean(y[zi == 0])
    else
      0
    prec1 <- 1 / s0^2 + n1 / s2_1_c
    prec2 <- 1 / s0^2 + n2 / s2_2_c
    mu1_c <- rnorm(1, (m0 / s0^2 + n1 * ybar1 / s2_1_c) / prec1, sqrt(1 / prec1))
    mu2_c <- rnorm(1, (m0 / s0^2 + n2 * ybar2 / s2_2_c) / prec2, sqrt(1 / prec2))
    
    # 3. Mixing proportion
    pi_c <- rbeta(1, a + n1, b + n2)
    
    # 4. Component-specific variances
    sse1 <- if (n1 > 0)
      sum((y[zi == 1] - mu1_c)^2)
    else
      0
    sse2 <- if (n2 > 0)
      sum((y[zi == 0] - mu2_c)^2)
    else
      0
    s2_1_c <- 1 / rgamma(1, shape = alpha0 + n1 / 2, rate = beta0 + sse1 / 2)
    s2_2_c <- 1 / rgamma(1, shape = alpha0 + n2 / 2, rate = beta0 + sse2 / 2)
    
    mu1[it]  <- mu1_c
    mu2[it]  <- mu2_c
    pi[it]   <- pi_c
    s2_1[it] <- s2_1_c
    s2_2[it] <- s2_2_c
    z1[it]   <- n1
  }
  
  out <- data.frame(
    mu1 = mu1,
    mu2 = mu2,
    pi = pi,
    s2_1 = s2_1,
    s2_2 = s2_2,
    n1 = z1
  )
  list(samples = out[seq.int(burn + 1, n_iter), , drop = FALSE],
       n = n,
       burn = burn)
}

chains_f_unequal <- lapply(seq_along(inits_f), function(i) {
  rgibbs_mix_unequal(
    y_faith,
    n_iter = 20000,
    burn = 5000,
    init = inits_f[[i]],
    seed = 261024400 + i
  )
})
names(chains_f_unequal) <- names(inits_f)

relabel_unequal <- function(ch) {
  s <- ch$samples
  flip <- s$mu1 > s$mu2
  data.frame(
    mu_lo  = ifelse(flip, s$mu2, s$mu1),
    mu_hi  = ifelse(flip, s$mu1, s$mu2),
    pi_lo  = ifelse(flip, 1 - s$pi, s$pi),
    s2_lo  = ifelse(flip, s$s2_2, s$s2_1),
    s2_hi  = ifelse(flip, s$s2_1, s$s2_2)
  )
}

rl_f_unequal <- do.call(rbind, lapply(chains_f_unequal, relabel_unequal))

tab_f_unequal <- data.frame(
  Parameter = c(
    "mu (short)",
    "mu (long)",
    "pi (short)",
    "sigma2 (short)",
    "sigma2 (long)"
  ),
  Median    = vapply(rl_f_unequal, median, numeric(1)),
  lo        = vapply(rl_f_unequal, quantile, numeric(1), probs = 0.025),
  hi        = vapply(rl_f_unequal, quantile, numeric(1), probs = 0.975),
  row.names = NULL
)
write.csv(
  tab_f_unequal,
  file.path(table_dir, "task2_faithful_unequal_variance_summary.csv"),
  row.names = FALSE
)

# Overlay Equal vs Unequal Variance Predictive Densities on Old Faithful
f_pred_unequal <- sapply(xs, function(x) {
  mean(
    rl_f_unequal$pi_lo * dnorm(x, rl_f_unequal$mu_lo, sqrt(rl_f_unequal$s2_lo)) + (1 - rl_f_unequal$pi_lo) * dnorm(x, rl_f_unequal$mu_hi, sqrt(rl_f_unequal$s2_hi))
  )
})

df_dens_comp <- rbind(
  data.frame(x = xs, density = f_pred, Model = "Equal variance (baseline)"),
  data.frame(x = xs, density = f_pred_unequal, Model = "Unequal variance (extended)")
)

fig_faithful_model_comparison <- ggplot(data.frame(y = y_faith), aes(y)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 30,
    fill = plot_palette[1],
    alpha = 0.25,
    colour = "white",
    linewidth = 0.2
  ) +
  geom_line(
    data = df_dens_comp,
    aes(
      x = x,
      y = density,
      colour = Model,
      linetype = Model
    ),
    inherit.aes = FALSE,
    linewidth = 0.65
  ) +
  scale_colour_manual(values = c(plot_palette[1], plot_palette[2])) +
  scale_linetype_manual(values = c(1, 2)) +
  labs(x = "Eruption duration (min)", y = "Density") +
  theme(legend.position = "top", legend.title = element_blank())

save_sci(
  fig_faithful_model_comparison,
  "task2_faithful_equal_vs_unequal_density",
  width = "double",
  height_mm = 75
)

## ---- 8. Save outputs and display the analysis summary ---------------------
# ---- 8.1 Record starting values, seeds and sampling settings ---------------
all_chains <- list(Separated = chains_sep,
                   Overlap = chains_ov,
                   Faithful = chains_f)
initial_values <- do.call(rbind, lapply(names(all_chains), function(nm)
  do.call(rbind, lapply(seq_along(all_chains[[nm]]), function(i) {
    ch <- all_chains[[nm]][[i]]
    data.frame(
      Setting = nm,
      Chain = names(all_chains[[nm]])[i],
      Seed = ch$seed,
      as.list(ch$init_used),
      Iterations = ch$n_iter,
      Burn = ch$burn,
      Sampler_seconds = ch$elapsed_seconds,
      row.names = NULL
    )
  }))))
write.csv(
  initial_values,
  file.path(table_dir, "task2_initial_values_and_seeds.csv"),
  row.names = FALSE
)
# ---- 8.2 Additional summaries for the complete results report --------------
chain_parameter_summary <- function(chains, setting) {
  do.call(rbind, lapply(names(chains), function(nm) {
    s <- relabel_chain(chains[[nm]])[, c("mu_lo", "mu_hi", "pi_lo", "s2")]
    data.frame(
      Setting = setting,
      Chain = nm,
      Parameter = names(s),
      Mean = vapply(s, mean, numeric(1)),
      Median = vapply(s, median, numeric(1)),
      SD = vapply(s, sd, numeric(1)),
      lo = vapply(s, quantile, numeric(1), probs = 0.025),
      hi = vapply(s, quantile, numeric(1), probs = 0.975),
      row.names = NULL
    )
  }))
}
chain_parameter_tab <- do.call(rbind, lapply(names(all_chains), function(nm)
  chain_parameter_summary(all_chains[[nm]], nm)))
write.csv(
  chain_parameter_tab,
  file.path(table_dir, "task2_chain_parameter_summary.csv"),
  row.names = FALSE
)

allocation_summary <- data.frame(
  Setting = c("Separated", "Overlap", "Faithful"),
  N = c(length(y_sep), length(y_ov), length(y_faith)),
  Ambiguous_count = c(
    sum(p_sep > 0.2 & p_sep < 0.8),
    sum(p_ov > 0.2 & p_ov < 0.8),
    sum(p_faith_alloc > 0.2 & p_faith_alloc < 0.8)
  ),
  Ambiguous_fraction = c(amb_sep, amb_ov, amb_faith),
  Mean_entropy = c(
    mean(alloc_sep$entropy),
    mean(alloc_ov$entropy),
    mean(alloc_f$entropy)
  ),
  MAP_accuracy = c(mean((p_sep >= .5) == z_true), mean((p_ov >= .5) == z_true), NA_real_),
  Brier_score = c(mean((p_sep - z_true)^2), mean((p_ov - z_true)^2), NA_real_)
)
write.csv(
  allocation_summary,
  file.path(table_dir, "task2_allocation_summary.csv"),
  row.names = FALSE
)

occupancy_tab <- do.call(rbind, lapply(names(all_chains), function(nm)
  do.call(rbind, lapply(names(all_chains[[nm]]), function(cn) {
    ch <- all_chains[[nm]][[cn]]
    s <- ch$samples
    data.frame(
      Setting = nm,
      Chain = cn,
      Empty_component_fraction = mean(s$n1 == 0 |
                                        s$n1 == ch$n),
      Extreme_weight_fraction = mean(s$pi < 0.01 |
                                       s$pi > 0.99)
    )
  }))))
write.csv(occupancy_tab,
          file.path(table_dir, "task2_component_occupancy.csv"),
          row.names = FALSE)

invariant_diagnostics <- function(chains, setting) {
  mats <- lapply(chains, function(ch) {
    s <- ch$samples
    cbind(
      mixture_mean = s$pi * s$mu1 + (1 - s$pi) * s$mu2,
      mixture_variance = s$s2 + s$pi * (1 - s$pi) * (s$mu1 - s$mu2)^2,
      separation = abs(s$mu1 - s$mu2)
    )
  })
  arr <- array(NA_real_,
               c(nrow(mats[[1]]), length(mats), 3),
               dimnames = list(NULL, names(chains), colnames(mats[[1]])))
  for (i in seq_along(mats))
    arr[, i, ] <- mats[[i]]
  d <- as.data.frame(
    posterior::summarise_draws(
      posterior::as_draws_array(arr),
      "mean",
      "median",
      "sd",
      "rhat",
      "ess_bulk",
      "ess_tail",
      "mcse_mean"
    )
  )
  names(d) <- c("Quantity",
                "Mean",
                "Median",
                "SD",
                "Rhat",
                "ESS",
                "ESS_tail",
                "MCSE")
  cbind(Setting = setting, d)
}
invariant_diag_tab <- do.call(rbind, lapply(names(all_chains), function(nm)
  invariant_diagnostics(all_chains[[nm]], nm)))
write.csv(
  invariant_diag_tab,
  file.path(table_dir, "task2_label_invariant_diagnostics.csv"),
  row.names = FALSE
)

# ---- 8.3 Save session information and posterior workspace ------------------
writeLines(capture.output(sessionInfo()),
           file.path(result_dir, "task2_sessionInfo.txt"))
save(
  chains_sep,
  chains_ov,
  chains_f,
  rl_sep,
  rl_ov,
  diag_tab,
  diag_raw_tab,
  chain_f,
  rl_f,
  tab_sep,
  tab_ov,
  tab_f,
  amb_sep,
  amb_ov,
  amb_faith,
  ess_f,
  diag_detailed,
  diag_raw_detailed,
  screen_tab,
  initial_values,
  label_tab,
  invariant_tab,
  chain_parameter_tab,
  allocation_summary,
  occupancy_tab,
  invariant_diag_tab,
  p_sep,
  p_ov,
  p_faith_alloc,
  alloc_sep,
  alloc_ov,
  alloc_f,
  density_sep,
  density_ov,
  y_sep,
  y_ov,
  y_faith,
  z_true,
  sw_sep,
  sw_ov,
  # Results from the team's Parts C/D additions.
  df_joint,
  tab_ov_relabel_critique,
  ppc_f,
  chains_f_unequal,
  rl_f_unequal,
  tab_f_unequal,
  xs,
  f_pred,
  f_pred_unequal,
  df_dens_comp,
  file = file.path(result_dir, "task2_sim.RData")
)
