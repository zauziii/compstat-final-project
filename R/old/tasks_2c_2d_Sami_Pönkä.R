# this is the file for additions to tasks 2c and 2d: 
# This section can be seamlessly added to the end of the main r code.

## ---- 8. Additions: Part C Relabeling Critique & Part D Faithful Follow-up --

# ---- 8.1 Part C: 2D Joint Posterior (mu1 vs mu2) & Empty-Component Effect ---
joint_trace_df <- function(chains, setting) {
  do.call(rbind, lapply(seq_along(chains), function(i) {
    ch <- chains[[i]]
    s  <- ch$samples
    keep <- unique(round(seq(1, nrow(s), length.out = min(2000, nrow(s)))))
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
  joint_trace_df(chains_ov,  "Overlapping components")
)
df_joint$setting <- factor(
  df_joint$setting,
  levels = c("Separated components", "Overlapping components")
)

fig_joint_posterior <- ggplot(df_joint, aes(x = mu1, y = mu2)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey40", linewidth = 0.4) +
  geom_point(aes(colour = chain, shape = empty), alpha = 0.55, size = 1.1) +
  facet_wrap(~ setting, scales = "free") +
  scale_colour_discrete(name = "Chain", limits = c("swap", "spread")) +
  scale_shape_manual(
    name   = "Empty component",
    values = c(`FALSE` = 16, `TRUE` = 4),
    labels = c(`FALSE` = "No (n1, n2 > 0)", `TRUE` = "Yes (prior draw)")
  ) +
  labs(
    x = expression(mu[1]),
    y = expression(mu[2])
  )

save_sci(
  fig_joint_posterior,
  "task2_joint_posterior_mu1_mu2",
  width = "double",
  height_mm = 80
)

# Compare all ordered draws vs. occupancy-conditioned ordered draws (min(n1, n2) >= 5)
relabel_filtered_summary <- function(chains, min_occ = 5) {
  ord_all <- do.call(rbind, lapply(chains, relabel_chain))
  n_obs   <- chains[[1]]$n
  keep    <- ord_all$n_lo >= min_occ & (n_obs - ord_all$n_lo) >= min_occ
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
    summarize_set(ord_sub, paste0("Occupancy-filtered (min n_k >= ", min_occ, ")"))
  )
}

tab_ov_relabel_critique <- relabel_filtered_summary(chains_ov, min_occ = 5)
write.csv(
  tab_ov_relabel_critique,
  file.path(table_dir, "task2_overlap_relabeling_critique.csv"),
  row.names = FALSE
)

# ---- 8.2 Part D: Formal Posterior Predictive Check (PPC) for Old Faithful ---
run_faithful_ppc <- function(rl_draws, y_obs, split_val = 3.0, n_rep = 2000, seed = 20261022) {
  set.seed(seed)
  n <- length(y_obs)
  idx <- sample(seq_len(nrow(rl_draws)), size = min(n_rep, nrow(rl_draws)), replace = FALSE)
  
  calc_stats <- function(y_vec) {
    s_short <- sd(y_vec[y_vec < split_val])
    s_long  <- sd(y_vec[y_vec >= split_val])
    c(sd_short = s_short, sd_long = s_long, sd_ratio = s_long / s_short)
  }
  
  obs_stats <- calc_stats(y_obs)
  rep_stats <- vapply(idx, function(m) {
    z_rep <- rbinom(n, 1, rl_draws$pi_lo[m])
    y_rep <- ifelse(
      z_rep == 1,
      rnorm(n, rl_draws$mu_lo[m], sqrt(rl_draws$s2[m])),
      rnorm(n, rl_draws$mu_hi[m], sqrt(rl_draws$s2[m]))
    )
    calc_stats(y_rep)
  }, numeric(3))
  
  rep_df <- as.data.frame(t(rep_stats))
  summary_df <- data.frame(
    Statistic = c("SD short (< 3 min)", "SD long (>= 3 min)", "Ratio (SD long / SD short)"),
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
  list(summary = summary_df, rep_draws = rep_df, obs = obs_stats)
}

ppc_f <- run_faithful_ppc(rl_f, y_faith)
write.csv(
  ppc_f$summary,
  file.path(table_dir, "task2_faithful_ppc_summary.csv"),
  row.names = FALSE
)

# ---- 8.3 Part D: Unequal-Variance Gibbs Sampler on Old Faithful -------------
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
  if (!is.null(seed)) set.seed(seed)
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
    ybar1 <- if (n1 > 0) mean(y[zi == 1]) else 0
    ybar2 <- if (n2 > 0) mean(y[zi == 0]) else 0
    prec1 <- 1 / s0^2 + n1 / s2_1_c
    prec2 <- 1 / s0^2 + n2 / s2_2_c
    mu1_c <- rnorm(1, (m0 / s0^2 + n1 * ybar1 / s2_1_c) / prec1, sqrt(1 / prec1))
    mu2_c <- rnorm(1, (m0 / s0^2 + n2 * ybar2 / s2_2_c) / prec2, sqrt(1 / prec2))
    
    # 3. Mixing proportion
    pi_c <- rbeta(1, a + n1, b + n2)
    
    # 4. Component-specific variances
    sse1 <- if (n1 > 0) sum((y[zi == 1] - mu1_c)^2) else 0
    sse2 <- if (n2 > 0) sum((y[zi == 0] - mu2_c)^2) else 0
    s2_1_c <- 1 / rgamma(1, shape = alpha0 + n1 / 2, rate = beta0 + sse1 / 2)
    s2_2_c <- 1 / rgamma(1, shape = alpha0 + n2 / 2, rate = beta0 + sse2 / 2)
    
    mu1[it]  <- mu1_c
    mu2[it]  <- mu2_c
    pi[it]   <- pi_c
    s2_1[it] <- s2_1_c
    s2_2[it] <- s2_2_c
    z1[it]   <- n1
  }
  
  out <- data.frame(mu1 = mu1, mu2 = mu2, pi = pi, s2_1 = s2_1, s2_2 = s2_2, n1 = z1)
  list(samples = out[seq.int(burn + 1, n_iter), , drop = FALSE], n = n, burn = burn)
}

chains_f_unequal <- lapply(seq_along(inits_f), function(i) {
  rgibbs_mix_unequal(y_faith, n_iter = 20000, burn = 5000, init = inits_f[[i]], seed = 400 + i)
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
  Parameter = c("mu (short)", "mu (long)", "pi (short)", "sigma2 (short)", "sigma2 (long)"),
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
  data.frame(x = xs, density = f_pred,         Model = "Equal variance (baseline)"),
  data.frame(x = xs, density = f_pred_unequal, Model = "Unequal variance (extended)")
)

fig_faithful_model_comparison <- ggplot(data.frame(y = y_faith), aes(y)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 30, fill = plot_palette[1], alpha = 0.25, colour = "white", linewidth = 0.2
  ) +
  geom_line(
    data = df_dens_comp,
    aes(x = x, y = density, colour = Model, linetype = Model),
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