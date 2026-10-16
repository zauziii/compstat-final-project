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
pkgs <- c("tidyverse", "scales", "ggpubr", "coda", "ggsci")

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

# --- publication ready style: global settings (fonts, sizes, colors) --------
theme_set(
  theme_pubr(base_size = 7, base_family = "Helvetica") +
    theme(
      palette.colour.discrete = "Okabe-Ito",
      palette.fill.discrete = "Okabe-Ito"
    )
)


## ---- 2. Part A: the Gibbs sampler -----------------------------------------
# Full-conditional Gibbs sampler for the 2-component Gaussian mixture.
rgibbs_mix <- function(y, n_iter = 5000, burn = 1000,
                       m0 = 0, s0 = 10, a = 1, b = 1,
                       alpha0 = 2, beta0 = 1,
                       init = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  y <- as.numeric(y)
  n <- length(y)
  mu1 <- mu2 <- pi <- s2 <- numeric(n_iter)
  z1 <- integer(n_iter)                       # n1 at each iteration
  
  if (is.null(init)) {                        # k-means start (10 restarts)
    km <- kmeans(y, centers = 2, nstart = 10)
    mu1_c <- km$centers[1]; mu2_c <- km$centers[2]
    pi_c  <- km$size[1] / n
    s2_c  <- mean((y - km$centers[km$cluster])^2)
  } else {
    mu1_c <- init$mu1
    mu2_c <- init$mu2
    pi_c  <- init$pi
    s2_c  <- init$s2
  }
  
  for (it in seq_len(n_iter)) {
    # 1. Allocations
    lik1 <- pi_c     * dnorm(y, mu1_c, sqrt(s2_c))
    lik2 <- (1-pi_c) * dnorm(y, mu2_c, sqrt(s2_c))
    p1   <- lik1 / (lik1 + lik2)
    zi   <- rbinom(n, 1, p1)
    n1   <- sum(zi == 1); n2 <- n - n1
    
    # 2. Component means
    ybar1 <- if (n1 > 0) mean(y[zi == 1]) else 0
    ybar2 <- if (n2 > 0) mean(y[zi == 0]) else 0
    prec1 <- 1/s0^2 + n1/s2_c
    mu1_c <- rnorm(1, (m0/s0^2 + n1*ybar1/s2_c)/prec1, sqrt(1/prec1))
    prec2 <- 1/s0^2 + n2/s2_c
    mu2_c <- rnorm(1, (m0/s0^2 + n2*ybar2/s2_c)/prec2, sqrt(1/prec2))
    
    # 3. Mixing proportion
    pi_c <- rbeta(1, a + n1, b + n2)
    
    # 4. Variance
    resid <- sum(ifelse(zi == 1, (y - mu1_c)^2, (y - mu2_c)^2))
    s2_c  <- 1 / rgamma(1, shape = alpha0 + n/2, rate = beta0 + resid/2)
    
    mu1[it] <- mu1_c; mu2[it] <- mu2_c
    pi[it]  <- pi_c;  s2[it]  <- s2_c
    z1[it]  <- n1
  }
  
  out <- data.frame(mu1 = mu1, mu2 = mu2, pi = pi, s2 = s2, n1 = z1)
  list(samples = out[(burn + 1):n_iter, , drop = FALSE],
       all = out, init_used = c(mu1 = mu1_c, mu2 = mu2_c, pi = pi_c, s2 = s2_c))
}


## ---- 4. Part B: simulation study ------------------------------------------
# Two settings: clearly separated vs substantially overlapping components
n <- 500
set.seed(20261016)
z_true <- rbinom(n, 1, 0.5)
y_sep <- ifelse(z_true == 1, rnorm(n, -2, 1),   rnorm(n, 2, 1))    # mu = +/- 2
y_ov  <- ifelse(z_true == 1, rnorm(n, -0.5, 1), rnorm(n, 0.5, 1))  # mu = +/- 0.5

# Three chains from meaningfully different starting values
inits <- list(
  swap   = list(mu1 = 2, mu2 = -2, pi = 0.9, s2 = 1),     # swapped labels
  spread = list(mu1 = 4, mu2 = -4, pi = 0.5, s2 = 4)      # spread means
)
run_chains <- function(y, inits, n_iter = 5000, burn = 1000) {
  lapply(seq_along(inits), function(i)
    rgibbs_mix(y, n_iter = n_iter, burn = burn, init = inits[[i]],
               seed = 100 + i))
}
chains_sep <- run_chains(y_sep, inits)
chains_ov  <- run_chains(y_ov,  inits)

# --- relabel by ordering the means (identifiability restriction) -----------
pool_relabelled <- function(chains) {
  do.call(rbind, lapply(chains, function(ch) {
    lo <- pmin(ch$samples$mu1, ch$samples$mu2)
    hi <- pmax(ch$samples$mu1, ch$samples$mu2)
    pi_lo <- ifelse(ch$samples$mu1 < ch$samples$mu2,
                    ch$samples$pi, 1 - ch$samples$pi)
    data.frame(mu_lo = lo, mu_hi = hi, pi_lo = pi_lo, s2 = ch$samples$s2)
  }))
}
rl_sep <- pool_relabelled(chains_sep)
rl_ov  <- pool_relabelled(chains_ov)

# --- posterior allocation probabilities (thinned) --------------------------
alloc_post <- function(ch, y, thin = 20) {
  s <- ch$samples
  keep <- seq(1, nrow(s), by = thin)
  P <- matrix(0, nrow = length(y), ncol = length(keep))
  for (j in seq_along(keep)) {
    i <- keep[j]
    l1 <- s$pi[i] * dnorm(y, s$mu1[i], sqrt(s$s2[i]))
    l2 <- (1 - s$pi[i]) * dnorm(y, s$mu2[i], sqrt(s$s2[i]))
    P[, j] <- l1 / (l1 + l2)
  }
  rowMeans(P)
}
p_sep <- alloc_post(chains_sep[[1]], y_sep)
p_ov  <- alloc_post(chains_ov[[1]], y_ov)
amb_sep <- mean(p_sep > 0.2 & p_sep < 0.8)
amb_ov  <- mean(p_ov  > 0.2 & p_ov  < 0.8)




