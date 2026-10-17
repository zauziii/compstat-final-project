# ============================================================================
#
#
#
#
# Log:
# 16.10.2026: folders are not being created in the correct path; add a 
#             helper function to locate the project root and set the 
#             working directory; finished the gibbs sampler; part B in progress.
# 17.10.2026: 
# ============================================================================


## ---- 0. project root ------------------------------------------------------
locate_root <- function() {
  # Get the script path from Rscript
  p <- grep("^--file=", commandArgs(), value = TRUE)
  p <- if (length(p)) sub("^--file=", "", p[1]) else ""
  
  # Prefer the active script in RStudio
  if (interactive() && requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    active <- rstudioapi::getSourceEditorContext()$path
    if (length(active) && nzchar(active)) p <- active
  }
  
  # Start from the script directory or working directory
  d <- if (nzchar(p) && file.exists(p)) {
    dirname(normalizePath(p, winslash = "/"))
  } else getwd()
  
  # Find the project root
  repeat {
    if (file.exists(file.path(d, "R", "task2_full_analysis.R"))) return(d)
    parent <- dirname(d)
    if (identical(parent, d)) stop("Could not find the project root.")
    d <- parent
  }
}
root <- locate_root()
setwd(root)


## ---- 0.5 output layout ----------------------------------------------------
data_dir <- "data"                  # raw inputs
table_dir <- "results/tables"       # tables ad csv
figure_dir <- "results/figures"     # figures as PDF + PNG
result_dir <- "results"             # RData workspace
report_dir <- "report"              # rendered report, quarto output
dir.create(table_dir, recursive = T, showWarnings = F)
dir.create(figure_dir, recursive = T, showWarnings = F)
dir.create(result_dir, recursive = T, showWarnings = F)
dir.create(report_dir, recursive = T, showWarnings = F)

## ---- 1. setup: libraries and visualization style --------------------------
# list all packages needed
pkgs <- c("tidyverse", "scales", "ggpubr", "coda", "ggsci", "posterior")

# get the names of all currently installed packages
installed <- rownames(installed.packages())
# find packages that are not installed yet
to_install <- setdiff(pkgs, installed)

# install only the missing packages
if (length(to_install) > 0) {
  install.packages(to_install)
}

# load the packages
# characer.only = TRUE lets library() accept a package name as a character string
invisible(lapply(pkgs, library, character.only = TRUE))

# publication ready style: global settings (fonts, sizes, colors)
theme_set(
  theme_pubr(base_size = 7, base_family = "Helvetica") +
    theme(
      palette.colour.discrete = "Okabe-Ito",
      palette.fill.discrete = "Okabe-Ito"
    )
)

# Export a figure at print size: 600 dpi PNG.
save_sci <- function(plot,
                     filename,
                     width = c("single", "double"),
                     height_mm = 68,
                     directory = figure_dir) {
  width <- match.arg(width)
  width_mm <- if (width == "single")
    89
  else
    183
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  
  out <- file.path(directory, paste0(filename, ".png"))
  ggsave(
    out,
    plot = plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    bg = "white",
    dpi = 600,
    device = grDevices::png
  )
  
  invisible(out)
}

## ---- 2. Part A: the Gibbs sampler -----------------------------------------
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


## ---- 4. Part B: simulation study ------------------------------------------
# Two settings: clearly separated vs substantially overlapping components
n <- 500
set.seed(20261016)
z_true <- rbinom(n, 1, 0.5)
y_sep <- ifelse(z_true == 1, rnorm(n, -2, 1), rnorm(n, 2, 1))    # mu = +/- 2
y_ov  <- ifelse(z_true == 1, rnorm(n, -0.5, 1), rnorm(n, 0.5, 1))  # mu = +/- 0.5

# Three chains from meaningfully different starting values
inits <- list(
  # k-means start
  km = NULL,
  # swapped labels
  swap   = list(
    mu1 = 2,
    mu2 = -2,
    pi = 0.9,
    s2 = 1
  ),
  # spread means
  spread = list(
    mu1 = 4,
    mu2 = -4,
    pi = 0.5,
    s2 = 4
  )
)
run_chains <- function(y,
                       inits,
                       n_iter = 5000,
                       burn = 1000,
                       seed_base = 100) {
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
                         seed_base = 100)
message("Overlapping simulation")
chains_ov  <- run_chains(y_ov,
                         inits,
                         n_iter = 60000,
                         burn = 15000,
                         seed_base = 200)

# --- relabel by ordering the means (identifiability restriction) -----------
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

# --- posterior allocation probabilities (thinned) --------------------------
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
    ambiguous_02_08 = p > .2 & p < .8,
    entropy = -safe * log(safe) - (1 - safe) * log1p(-safe)
  )
}
p_sep <- alloc_chains(chains_sep, y_sep)
p_ov  <- alloc_chains(chains_ov, y_ov)
alloc_sep <- allocation_table(y_sep, p_sep)
alloc_ov <- allocation_table(y_ov, p_ov)
write.csv(alloc_sep,
          file.path(table_dir, "allocations_separated.csv"),
          row.names = FALSE)
write.csv(alloc_ov,
          file.path(table_dir, "allocations_overlap.csv"),
          row.names = FALSE)
amb_sep <- mean(p_sep > 0.2 & p_sep < 0.8)
amb_ov  <- mean(p_ov  > 0.2 & p_ov  < 0.8)

# --- Table 1: recovery of the data-generating mixture -----------------------
mk_rec <- function(truth, mu_lo, mu_hi, pi_lo, s2) {
  data.frame(
    Parameter = c("mu1 (lower)", "mu2 (upper)", "pi", "sigma2"),
    Truth = truth,
    Median = c(median(mu_lo), median(mu_hi), median(pi_lo), median(s2)),
    lo = c(
      quantile(mu_lo, .025),
      quantile(mu_hi, .025),
      quantile(pi_lo, .025),
      quantile(s2, .025)
    ),
    hi = c(
      quantile(mu_lo, .975),
      quantile(mu_hi, .975),
      quantile(pi_lo, .975),
      quantile(s2, .975)
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
write.csv(tab_sep,
          file.path(table_dir, "table1_recovery_separated.csv"),
          row.names = FALSE)
write.csv(tab_ov,
          file.path(table_dir, "table1_recovery_overlap.csv"),
          row.names = FALSE)

## ---- 5. Figure 1: trace plots ---------------------------------------------
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
  geom_line(linewidth = 0.25, alpha = 0.75) +
  scale_color_manual(
    name = NULL,
    values = c(mu1 = "#D55E00", mu2 = "#0072B2"),
    labels = c(mu1 = expression(mu[1]), mu2 = expression(mu[2]))
  ) +
  scale_linetype_discrete(name = "Chain") +
  labs(x = NULL, y = "MCMC draw") +
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
  geom_line(linewidth = 0.25, alpha = 0.75) +
  scale_color_manual(
    name = NULL,
    values = c(mu1 = "#D55E00", mu2 = "#0072B2"),
    labels = c(mu1 = expression(mu[1]), mu2 = expression(mu[2]))
  ) +
  scale_linetype_discrete(name = "Chain") +
  labs(x = "Iteration (post burn-in)", y = "MCMC draw") +
  theme(legend.position = "none")

fig1 <- ggarrange(
  p_trace_sep,
  p_trace_ov,
  ncol = 1,
  labels = c("a", "b"),
  font.label = list(size = 8, face = "bold"),
  common.legend = TRUE,
  legend = "right",
  align = "v"
)

## ---- 6. Table 2: diagnostics (ESS, Gelman-Rubin) -------------------------
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
write.csv(diag_tab, file.path(table_dir, "table2_diagnostics.csv"))
write.csv(diag_raw_tab,
          file.path(table_dir, "table2_diagnostics_raw.csv"))
write.csv(
  diag_detailed,
  file.path(table_dir, "table2_diagnostics_detailed.csv"),
  row.names = FALSE
)
write.csv(
  diag_raw_detailed,
  file.path(table_dir, "table2_diagnostics_raw_detailed.csv"),
  row.names = FALSE
)

diagnostic_screen <- function(d, setting) {
  good <- all(is.finite(d$Rhat)) && max(d$Rhat) <= 1.01 &&
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

## ---- 7. Part C: label switching -------------------------------------------

