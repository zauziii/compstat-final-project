# ============================================================================
# MAST32001 Computational Statistics - Task 1
# Ladle estimator for the dimension of the signal subspace (SOBI / AMUSE).
#
# Log:
# 20.10.2026: create and finish this file.
# 21.10.2026: move the figures out of the run functions.
# 22.10.2026: save the block-length check and update older RData files.
# ============================================================================

## ---- 0. shared project setup ----------------------------------------------
# Run from the project root; the shared helper supplies paths, packages,
# the global theme, palette and PNG export function.
source(file.path("R", "project_setup.R"))
ensure_packages(c("JADE"))

## ---- 1. Part A: whitening and autocovariance matrices ---------------------
# Conventions used throughout:
#   * X is a T x p matrix, rows are time points, columns are series.
#   * Whitening: X_st = (X - 1_T * mu') * S^(-1/2), where mu is the column mean
#     and S = (1/T) * (X - mu)'(X - mu) is the covariance matrix scaled by 1/T.
#     S^(-1/2) = V D^(-1/2) V' from the eigendecomposition S = V D V'.
#   * Autocovariance at lag tau>0:  A_tau = (1/(T - tau)) * sum_{t=1}^{T-tau}
#     x_st(t) x_st(t+tau)', where x_st(t) is row t of the whitened series.
#     Both methods use its symmetric part: (A_tau + A_tau') / 2.
whiten <- function(X) {
  X <- as.matrix(X)
  n <- nrow(X)
  if (!is.numeric(X) ||
      n < 2L || ncol(X) < 1L || any(!is.finite(X)))
    stop("X must be a finite numeric matrix with at least two rows")
  mu <- colMeans(X)
  Xc <- sweep(X, 2L, mu)
  S <- crossprod(Xc) / n
  ev <- eigen(S, symmetric = TRUE)
  tol <- .Machine$double.eps * max(n, ncol(X)) * max(ev$values)
  if (max(ev$values) <= 0 || any(ev$values <= tol)) {
    stop(structure(
      list(message = "Sample covariance is singular or numerically singular", call = NULL),
      class = c("task1_singular_covariance", "error", "condition")
    ))
  }
  half_inv <- ev$vectors %*% diag(1 / sqrt(ev$values), nrow = ncol(X)) %*% t(ev$vectors)
  Xc %*% half_inv
}

autocov <- function(Xst, tau = 1L) {
  Xst <- as.matrix(Xst)
  n <- nrow(Xst)
  if (length(tau) != 1L || !is.finite(tau) || tau != floor(tau) ||
      tau < 1L || tau >= n)
    stop("lag tau must satisfy 1 <= tau < n")
  if (!is.numeric(Xst) ||
      any(!is.finite(Xst)))
    stop("Xst must be finite and numeric")
  A <- crossprod(Xst[seq_len(n - tau), , drop = FALSE], Xst[seq.int(tau + 1L, n), , drop = FALSE]) / (n - tau)
  (A + t(A)) / 2
}

# Numerically verify that whitening makes the sample covariance the identity.
whiten_check <- function(X) {
  Xst <- whiten(X)
  S <- crossprod(sweep(Xst, 2L, colMeans(Xst))) / nrow(Xst)
  off <- S - diag(ncol(S))
  data.frame(
    max_abs_diag_deviation = max(abs(diag(S) - 1)),
    max_abs_off_diagonal = if (ncol(S) == 1L)
      0
    else
      max(abs(off[row(off) != col(off)])),
    frobenius_error = sqrt(sum(off^2))
  )
}

## ---- 2. Part B: bootstrap for dependent observations ----------------------
# Fixed-block bootstrap (Kuensch 1989): use overlapping blocks of fixed
# length B; sample h = ceiling(T/B) block start indices with replacement from
# {1, ..., T - B + 1} and concatenate the corresponding blocks, truncating to
# the first T rows. Tuning parameter: the block length B.
fixed_block_bootstrap <- function(X, B) {
  X <- as.matrix(X)
  n <- nrow(X)
  if (length(B) != 1L || !is.finite(B) || B != floor(B) ||
      B < 1L || B > n)
    stop("block length B must satisfy 1 <= B <= T")
  starts <- sample.int(n - B + 1L, size = ceiling(n / B), replace = TRUE)
  idx <- unlist(lapply(starts, function(s)
    seq.int(s, length.out = B)))
  X[idx[seq_len(n)], , drop = FALSE]
}

# Stationary bootstrap (Politis & Romano 1994): each block has a random length
# drawn from a geometric distribution with success probability p, so the
# expected block length is 1/p. The time series is treated as a circle: a
# block starting near the end wraps around to the beginning, as in the
# original construction. Block starts are drawn uniformly from the observed
# time points and blocks are concatenated until the replicate is at least T
# long, then truncated. Tuning parameter: p (or the expected length).
stationary_bootstrap <- function(X, p) {
  X <- as.matrix(X)
  n <- nrow(X)
  if (n < 1L || length(p) != 1L || !is.finite(p) || p <= 0 || p > 1)
    stop("probability p must lie in (0, 1], with at least one observed row")
  idx <- integer(0)
  while (length(idx) < n) {
    start <- sample.int(n, size = 1L)
    len <- stats::rgeom(1, prob = p) + 1L       # untruncated length has E(len) = 1/p
    len <- min(len, n - length(idx))          # discard only the unused final tail
    idx <- c(idx, ((start - 1L + seq_len(len) - 1L) %% n) + 1L)
  }
  X[idx[seq_len(n)], , drop = FALSE]
}

# How well does each resampling scheme preserve the autocovariance of the
# original series? Used for the graphical comparison in the report.
bootstrap_demo <- function(x,
                           B = 20,
                           p = 1 / 20,
                           nboot = 200,
                           max_lag = 40) {
  if (!is.numeric(x) ||
      any(!is.finite(x)) || length(x) < 3L || stats::sd(x) <= 0)
    stop("x must be a nonconstant finite numeric series")
  if (length(nboot) != 1L ||
      !is.finite(nboot) || nboot < 2 || nboot != floor(nboot))
    stop("nboot must be an integer of at least two")
  if (length(max_lag) != 1L || !is.finite(max_lag) || max_lag < 1 ||
      max_lag >= length(x) || max_lag != floor(max_lag))
    stop("max_lag must be an integer between 1 and length(x)-1")
  acf_orig <- stats::acf(x, lag.max = max_lag, plot = FALSE)$acf[-1L]
  grab <- function(fun) {
    ac <- matrix(NA_real_, nboot, max_lag)
    for (j in seq_len(nboot)) {
      xb <- fun()
      a <- stats::acf(xb, lag.max = max_lag, plot = FALSE)$acf[-1L]
      ac[j, ] <- if (length(a) == max_lag)
        a
      else
        c(a, rep(NA_real_, max_lag - length(a)))
    }
    colMeans(ac, na.rm = TRUE)
  }
  X <- matrix(x, ncol = 1L)
  id <- grab(function()
    X[sample.int(length(x), replace = TRUE), , drop = FALSE])
  fb <- grab(function()
    fixed_block_bootstrap(X, B)[, 1L])
  sb <- grab(function()
    stationary_bootstrap(X, p)[, 1L])
  data.frame(
    lag = seq_len(max_lag),
    original = acf_orig,
    iid = id,
    fixed_block = fb,
    stationary = sb
  )
}

## ---- 3. Part C: the ladle quantities --------------------------------------
# Returns the p x p orthogonal joint diagonalizer U for the autocovariance
# matrices at the given lags. AMUSE uses a single lag and an ordinary
# eigendecomposition; SOBI uses several lags and JADE::rjd() (the assignment
# explicitly allows rjd for more than two matrices). If rjd fails to converge
# on a particular sample, NULL is returned. Singular bootstrap covariances
# also count as numerical failures. Other input errors are not suppressed.
# Crucially, U's columns are sorted by DECREASING SQUARED lag strength.
joint_diagonalizer <- function(X,
                               method = c("SOBI", "AMUSE"),
                               lags = 1:12,
                               tau = 1L) {
  method <- match.arg(method)
  Xst <- tryCatch(
    whiten(X),
    task1_singular_covariance = function(e)
      NULL
  )
  if (is.null(Xst))
    return(NULL)
  p <- ncol(Xst)
  lags_use <- if (method == "AMUSE")
    tau
  else
    lags
  if (!length(lags_use) ||
      any(!is.finite(lags_use)) || anyDuplicated(lags_use))
    stop("lags must be a nonempty vector of distinct finite integers")
  mats <- array(NA_real_, dim = c(p, p, length(lags_use)))
  for (j in seq_along(lags_use))
    mats[, , j] <- autocov(Xst, lags_use[j])
  if (method == "AMUSE") {
    U <- eigen(mats[, , 1L], symmetric = TRUE)$vectors
  } else {
    fit <- tryCatch(
      JADE::rjd(mats, maxiter = 1000),
      error = function(e)
        NULL
    )
    if (is.null(fit))
      return(NULL)
    U <- fit$V
  }
  strength <- Reduce(`+`, lapply(seq_along(lags_use), function(j)
    diag(crossprod(U, mats[, , j] %*% U))^2))
  U[, order(strength, decreasing = TRUE), drop = FALSE]
}

# Squared joint diagonalization eigenvalues, ordered decreasingly: the
# diagonal elements of Lambda = sum_tau diag(U' A_tau U)^2 (SOBI) or
# Lambda = diag(U' A_tau U)^2 (AMUSE). sort_values=FALSE retains the
# supplied column order, as required when pairing these scores with U.
squared_eigenvalues <- function(X,
                                U,
                                method = c("SOBI", "AMUSE"),
                                lags = 1:12,
                                tau = 1L,
                                sort_values = TRUE) {
  method <- match.arg(method)
  Xst <- whiten(X)
  lags_use <- if (method == "AMUSE")
    tau
  else
    lags
  lam2 <- Reduce(`+`, lapply(lags_use, function(lag) {
    A <- autocov(Xst, lag)
    diag(t(U) %*% A %*% U)^2
  }))
  if (sort_values)
    sort(lam2, decreasing = TRUE)
  else
    lam2
}

# The ladle criterion of Nordhausen & Virta (2018):
#   f_0(k)   = (1/n_boot) sum_j {1 - |det(Bhat_k' B*_{k,j})|},  f_0(0) = 0
#   f(k)     = f_0(k) / (1 + sum_{l=0}^{p-1} f_0(l))
#   phi(k)   = lambda_{k+1}^2 / (1 + sum_{j=1}^{p} lambda_j^2), k=0,...,p-1
#   g(k)     = f(k) + phi(k),     qhat = argmin_{k=0..p-1} g(k)
# Bhat_k and B*_{k,j} are the first k columns of the joint diagonalizer
# computed on the original sample and on the j-th bootstrap sample. Bootstrap
# samples that fail numerically are recorded. At least min_boot_fraction
# (default 90%) must succeed; otherwise qhat is NA and status explains why.
# lambda2 already contains squared strengths: it must not be squared again.
ladle_estimate <- function(X,
                           method = c("SOBI", "AMUSE"),
                           boot = c("F", "S"),
                           block = NULL,
                           p_boot = NULL,
                           nboot = 200,
                           lags = 1:12,
                           tau = 1L,
                           seed = NULL,
                           min_boot_fraction = 0.9) {
  method <- match.arg(method)
  boot <- match.arg(boot)
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  if (p < 2L || !is.numeric(X) || any(!is.finite(X)))
    stop("X must contain at least two finite numeric series")
  if (length(nboot) != 1L ||
      !is.finite(nboot) || nboot < 2L || nboot != floor(nboot))
    stop("nboot must be an integer of at least two")
  if (length(min_boot_fraction) != 1L ||
      !is.finite(min_boot_fraction) ||
      min_boot_fraction <= 0 || min_boot_fraction > 1)
    stop("min_boot_fraction must lie in (0, 1]")
  if (boot == "F" && (length(block) != 1L || !is.finite(block) ||
                      block != floor(block) ||
                      block < 1L || block > n))
    stop("For F, block must be an integer between 1 and nrow(X)")
  if (boot == "S" && (length(p_boot) != 1L || !is.finite(p_boot) ||
                      p_boot <= 0 || p_boot > 1))
    stop("For S, p_boot must lie in (0, 1]")
  if (!is.null(seed))
    set.seed(seed)
  # Return the same fields on success and failure so no failed fit disappears.
  out <- list(
    k = 0:(p - 1L),
    f = rep(NA_real_, p),
    phi = rep(NA_real_, p),
    g = rep(NA_real_, p),
    qhat = NA_integer_,
    lambda2 = rep(NA_real_, p),
    dets = matrix(NA_real_, nboot, p - 1L),
    Bhat = NULL,
    n_boot_used = 0L,
    n_boot_attempted = 0L,
    n_boot_failed = 0L,
    method = method,
    boot = boot,
    block = block,
    p_boot = p_boot,
    nboot = nboot,
    lags = lags,
    tau = tau,
    seed = seed,
    min_boot_fraction = min_boot_fraction,
    g_gap = NA_real_,
    status = "original_fit_failed"
  )
  # Original-sample joint diagonalizer and squared eigenvalues.
  Bhat <- joint_diagonalizer(X, method, lags, tau)
  if (is.null(Bhat))
    return(out)
  lam2 <- squared_eigenvalues(X, Bhat, method, lags, tau, sort_values = FALSE)
  out$Bhat <- Bhat
  out$lambda2 <- lam2
  out$phi <- lam2 / (1 + sum(lam2))
  # Bootstrap replicates of the diagonalizer; k runs over 1, ..., p-1.
  dets <- matrix(NA_real_, nrow = nboot, ncol = p - 1L)
  n_ok <- 0L
  for (j in seq_len(nboot)) {
    Xb <- if (boot == "F")
      fixed_block_bootstrap(X, block)
    else
      stationary_bootstrap(X, p_boot)
    Bstar <- joint_diagonalizer(Xb, method, lags, tau)
    if (is.null(Bstar))
      next
    n_ok <- n_ok + 1L
    for (k in seq_len(p - 1L)) {
      A <- crossprod(Bhat[, seq_len(k), drop = FALSE], Bstar[, seq_len(k), drop = FALSE])
      dets[j, k] <- min(1, max(0, abs(det(A))))  # constrain roundoff to [0, 1]
    }
  }
  out$dets <- dets
  out$n_boot_used <- n_ok
  out$n_boot_attempted <- nboot
  out$n_boot_failed <- nboot - n_ok
  if (n_ok < max(2L, ceiling(nboot * min_boot_fraction))) {
    out$status <- "insufficient_bootstrap_success"
    return(out)
  }
  f0 <- c(0, colMeans(1 - dets, na.rm = TRUE))
  f <- f0 / (1 + sum(f0))
  phi <- lam2 / (1 + sum(lam2))
  g <- f + phi
  out$f <- f
  out$phi <- phi
  out$g <- g
  out$qhat <- which.min(g) - 1L
  out$g_gap <- diff(sort(g)[1:2])
  out$status <- if (n_ok == nboot)
    "ok"
  else
    "ok_with_bootstrap_failures"
  out
}

## ---- 4. Part D: simulation study (Setting 1 of the paper) -------------------
# Setting 1: p = 10, q = 5; the five signal components are
#   ARMA(0.6; (0,-0.4)), ARMA((0.4,0.1,0.3); (0.2,0.4)), AR(0.4),
#   MA(0.5,0.3,0.1,0.2), AR(0.1,0.1,0.3)
# with Gaussian innovations; the five noise components are Gaussian white
# noise. The mixing is the identity Omega = I_10 (affine equivariance).
# Innovations have variance one; marginal component scales are removed by
# whitening, equivalently absorbing them into a diagonal mixing matrix.
# A burn-in of n.start = 1000 observations removes the transient of arima.sim
# so that each generated signal is (close to) its stationary regime.
signal_models <- list(
  list(ar = c(0.6), ma = c(0, -0.4)),
  list(ar = c(0.4, 0.1, 0.3), ma = c(0.2, 0.4)),
  list(ar = c(0.4)),
  list(ar = numeric(0), ma = c(0.5, 0.3, 0.1, 0.2)),
  list(ar = c(0.1, 0.1, 0.3))
)

generate_setting1 <- function(T, seed = NULL) {
  if (length(T) != 1L || !is.finite(T) || T < 2L || T != floor(T))
    stop("T must be an integer of at least two")
  if (!is.null(seed))
    set.seed(seed)
  sig <- vapply(signal_models, function(m) {
    as.numeric(stats::arima.sim(
      model = m,
      n = T,
      n.start = 1000L,
      # burn-in: discard the transient
      innov = stats::rnorm(T)
    ))
  }, numeric(T))
  noise <- matrix(stats::rnorm(T * 5), T, 5)
  cbind(sig, noise)
}

# One Monte Carlo repetition: estimate the dimension for one (method, boot,
# block length, sample size) combination. NA marks any unusable estimate.
# details=TRUE also returns failure status and bootstrap success counts.
mc_once <- function(T,
                    method,
                    boot,
                    block,
                    nboot = 200,
                    seed = NULL,
                    boot_seed = NULL,
                    details = FALSE) {
  X <- generate_setting1(T, seed = seed)
  fit <- if (boot == "F")
    ladle_estimate(X,
                   method,
                   "F",
                   block = block,
                   nboot = nboot,
                   seed = boot_seed)
  else
    ladle_estimate(
      X,
      method,
      "S",
      p_boot = 1 / block,
      nboot = nboot,
      seed = boot_seed
    )
  if (!details)
    return(fit$qhat)
  data.frame(
    qhat = fit$qhat,
    status = fit$status,
    n_boot_used = fit$n_boot_used,
    n_boot_attempted = fit$n_boot_attempted,
    n_boot_failed = fit$n_boot_failed,
    nboot = fit$nboot,
    data_seed = if (is.null(seed))
      NA_integer_
    else
      seed,
    bootstrap_seed = if (is.null(boot_seed))
      NA_integer_
    else
      boot_seed
  )
}

# Keep all attempted repetitions, including cells in which every fit failed.
# Accuracy uses all attempts; the conditional success rate is reported
# separately. Wilson intervals remain meaningful when accuracy is 0 or 1.
summarise_simulation <- function(sim_results, q_true = 5L) {
  keys <- c("T", "method", "boot", "block")
  groups <- split(sim_results, interaction(sim_results[keys], drop = TRUE))
  do.call(rbind, lapply(groups, function(d) {
    q <- d$qhat
    valid <- is.finite(q)
    n <- length(q)
    n_valid <- sum(valid)
    n_correct <- sum(q[valid] == q_true)
    prob <- n_correct / n
    z <- stats::qnorm(0.975)
    centre <- (prob + z^2 / (2 * n)) / (1 + z^2 / n)
    half <- z * sqrt(prob * (1 - prob) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
    cbind(
      d[1L, keys, drop = FALSE],
      data.frame(
        n_attempted = n,
        n_valid = n_valid,
        n_failed = n - n_valid,
        failure_rate = 1 - n_valid / n,
        mean_qhat = if (n_valid)
          mean(q[valid])
        else
          NA_real_,
        p_qhat5 = prob,
        sd_qhat = if (n_valid >= 2L)
          stats::sd(q[valid])
        else
          NA_real_,
        p_qhat5_given_valid = if (n_valid)
          n_correct / n_valid
        else
          NA_real_,
        mcse = sqrt(prob * (1 - prob) / n),
        ci_low = max(0, centre - half),
        ci_high = min(1, centre + half),
        bias_valid = if (n_valid)
          mean(q[valid] - q_true)
        else
          NA_real_,
        rmse_valid = if (n_valid)
          sqrt(mean((q[valid] - q_true)^2))
        else
          NA_real_,
        p_under = sum(q[valid] < q_true) / n,
        p_over = sum(q[valid] > q_true) / n,
        n_boot_used_total = sum(d$n_boot_used, na.rm = TRUE),
        n_boot_attempted_total = sum(d$n_boot_attempted, na.rm = TRUE),
        n_boot_failed_total = sum(d$n_boot_failed, na.rm = TRUE),
        n_boot_counts_unknown = sum(is.na(d$n_boot_attempted))
      )
    )
  }))
}

# Paired comparison of SOBI and AMUSE on the same (T, repetition) datasets:
# the difference of the correct-dimension indicators, averaged over the
# datasets shared by both methods, with the paired MCSE of that difference.
# Requires sim_results with the repetition column (produced by run_task1).
paired_method_comparison <- function(sim_results, q_true = 5L) {
  keys <- c("T", "repetition", "boot", "block")
  sobi <- sim_results[sim_results$method == "SOBI", c(keys, "qhat")]
  amuse <- sim_results[sim_results$method == "AMUSE", c(keys, "qhat")]
  merged <- merge(sobi,
                  amuse,
                  by = keys,
                  suffixes = c("_SOBI", "_AMUSE"))
  ok <- function(q)
    as.integer(!is.na(q) & q == q_true)
  merged$success_SOBI <- ok(merged$qhat_SOBI)
  merged$success_AMUSE <- ok(merged$qhat_AMUSE)
  merged$difference <- merged$success_SOBI - merged$success_AMUSE
  groups <- split(merged, interaction(merged[c("T", "boot", "block")], drop = TRUE))
  do.call(rbind, lapply(groups, function(d) {
    data.frame(
      T = d$T[1L],
      boot = d$boot[1L],
      block = d$block[1L],
      paired_datasets = nrow(d),
      SOBI_minus_AMUSE_correct_rate = mean(d$difference),
      paired_mcse = if (nrow(d) > 1L)
        stats::sd(d$difference) / sqrt(nrow(d))
      else
        NA_real_
    )
  }))
}

## ---- 5. Part E: application to the course data ------------------------------
read_timeseries <- function(path) {
  if (!file.exists(path))
    stop("Course data file not found: ", path)
  X <- as.matrix(utils::read.table(path, header = FALSE))
  if (!is.numeric(X) ||
      nrow(X) < 3L || ncol(X) < 2L || any(!is.finite(X)))
    stop("Course data must contain at least two finite numeric series")
  colnames(X) <- paste0("V", seq_len(ncol(X)))
  X
}

## ---- 6. analysis and output -------------------------------------------------
run_task1 <- function(root = locate_root()) {
  if (is.null(root) ||
      !dir.exists(root))
    stop("Could not find the project root")
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(root)
  for (d in c(table_dir, figure_dir, result_dir))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  write_table <- function(x, name)
    utils::write.csv(x, file.path(table_dir, paste0(name, ".csv")), row.names = FALSE)
  # The defaults reproduce the original design. Environment overrides allow a
  # small validation run without editing the analysis or labelling it full MC.
  get_integer <- function(name, default, minimum = 1L) {
    text <- Sys.getenv(name, unset = as.character(default))
    value <- suppressWarnings(as.numeric(text))
    if (length(value) != 1L ||
        !is.finite(value) || value < minimum ||
        value != floor(value) || value > .Machine$integer.max)
      stop(name, " must be an integer >= ", minimum)
    as.integer(value)
  }
  get_grid <- function(name, default, minimum) {
    text <- Sys.getenv(name, unset = paste(default, collapse = ","))
    value <- suppressWarnings(as.numeric(strsplit(text, ",", fixed = TRUE)[[1L]]))
    if (!length(value) ||
        any(!is.finite(value)) || any(value < minimum) ||
        any(value != floor(value)) || anyDuplicated(value))
      stop(name,
           " must contain distinct comma-separated integers >= ",
           minimum)
    value
  }
  nboot_demo <- get_integer("TASK1_NBOOT_DEMO", 200L, 2L)
  nboot_example <- get_integer("TASK1_NBOOT_EXAMPLE", 200L, 2L)
  nboot_app <- get_integer("TASK1_NBOOT_APP", 200L, 2L)
  blocks <- get_grid("TASK1_BLOCKS", c(10, 20, 40), 1L)
  Ts <- get_grid("TASK1_TS", c(200, 500, 1000, 2000), 13L)
  if (any(blocks > min(c(Ts, 500L))))
    stop("Block lengths must not exceed the shortest series")
  nrep <- get_integer("TASK1_NREP", 200L, 2L)
  nboot_screen <- get_integer("TASK1_NBOOT", 100L, 2L)
  detected <- parallel::detectCores()
  default_cores <- if (is.na(detected))
    1L
  else
    max(1L, detected - 1L)
  ncores <- min(nrep, get_integer("TASK1_NCORES", default_cores))
  if (.Platform$OS.type == "windows" && ncores > 1L) {
    message("Windows: using one worker because mclapply does not support forking here.")
    ncores <- 1L
  }
  run_settings <- list(
    nrep = nrep,
    nboot_screen = nboot_screen,
    nboot_demo = nboot_demo,
    nboot_example = nboot_example,
    nboot_app = nboot_app,
    blocks = blocks,
    Ts = Ts,
    ncores = ncores,
    min_boot_fraction = 0.9,
    autocovariance = "symmetric; denominator T-lag"
  )
  write_table(data.frame(
    setting = names(run_settings),
    value = vapply(run_settings, function(x)
      paste(x, collapse = ","), "")
  ),
  "task1_run_settings")
  
  ## ---- Part A: whitening check ---------------------------------------------
  Xdemo <- generate_setting1(500, seed = 20261016)
  whiten_num <- whiten_check(Xdemo)
  message("Part A: whitening check on Setting-1 data (T=500)")
  print(whiten_num, row.names = FALSE)
  write_table(whiten_num, "task1_whitening_check")
  
  ## ---- Part B: bootstrap demonstration --------------------------------------
  x1 <- generate_setting1(500, seed = 20261017)[, 1L]  # ARMA(0.6; (0,-0.4)) signal
  set.seed(20261017)
  boot_demo <- do.call(rbind, lapply(blocks, function(b) {
    cbind(block = b,
          bootstrap_demo(
            x1,
            B = b,
            p = 1 / b,
            nboot = nboot_demo,
            max_lag = 40
          ))
  }))
  write_table(boot_demo, "task1_bootstrap_demo")
  
  ## ---- Part C: a single ladle plot (Setting 1, SOBI, stationary, b=20) -----
  set.seed(20261018)
  Xc <- generate_setting1(1000)
  ladle_example <- ladle_estimate(
    Xc,
    method = "SOBI",
    boot = "S",
    p_boot = 1 / 20,
    nboot = nboot_example,
    seed = 20261018
  )
  message(
    "Part C example (SOBI, S, E[b]=20, T=1000): qhat = ",
    ladle_example$qhat,
    " (",
    ladle_example$status,
    "; ",
    ladle_example$n_boot_used,
    "/",
    nboot_example,
    " bootstrap fits used)"
  )
  write_table(
    data.frame(
      k = ladle_example$k,
      f = ladle_example$f,
      phi = ladle_example$phi,
      g = ladle_example$g
    ),
    "task1_ladle_example_criteria"
  )
  
  ## ---- Part D: Monte Carlo study -------------------------------------------
  # Grid: method x bootstrap x block length x sample size. Repetitions and
  # nboot are chosen so that the study supports the qualitative conclusions
  # while the computational effort stays sensible (the paper uses 2000
  # repetitions and 200 bootstrap samples; here we use 200 and 100 and state
  # this clearly in the report). Use forked workers where supported.
  # The same dataset is used for all configurations within a (T, repetition)
  # pair; bootstrap seeds are separate and all seeds are saved in the CSV.
  methods <- c("SOBI", "AMUSE")
  boots <- c("F", "S")
  grid <- expand.grid(
    T = Ts,
    method = methods,
    boot = boots,
    block = blocks,
    stringsAsFactors = FALSE
  )
  message(
    "Part D: Monte Carlo over ",
    nrow(grid),
    " cells x ",
    nrep,
    " repetitions on ",
    ncores,
    " cores (this takes a while)"
  )
  sim_list <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    cell <- grid[i, ]
    seeds <- 100000 * match(cell$T, Ts) + seq_len(nrep)
    boot_seeds <- 10000000 + 100000 * i + seq_len(nrep)
    if (max(c(seeds, boot_seeds)) > .Machine$integer.max)
      stop("Seed range exceeds R's integer limit")
    reps <- parallel::mclapply(seq_len(nrep), function(r) {
      tryCatch(
        mc_once(
          cell$T,
          cell$method,
          cell$boot,
          cell$block,
          nboot = nboot_screen,
          seed = seeds[r],
          boot_seed = boot_seeds[r],
          details = TRUE
        ),
        error = function(e)
          data.frame(
            qhat = NA_integer_,
            status = paste("error:", conditionMessage(e)),
            n_boot_used = NA_integer_,
            n_boot_attempted = NA_integer_,
            n_boot_failed = NA_integer_,
            nboot = nboot_screen,
            data_seed = seeds[r],
            bootstrap_seed = boot_seeds[r]
          )
      )
    }, mc.cores = ncores, mc.set.seed = FALSE)
    if (any(vapply(reps, function(x)
      inherits(x, "try-error") || !is.data.frame(x), logical(1))))
      stop(
        "A parallel worker failed before returning its diagnostics; rerun with TASK1_NCORES=1"
      )
    draws <- do.call(rbind, reps)
    qs <- draws$qhat
    sim_list[[i]] <- cbind(
      data.frame(
        T = cell$T,
        method = as.character(cell$method),
        boot = as.character(cell$boot),
        block = cell$block,
        repetition = seq_len(nrep)
      ),
      draws
    )
    # Update the raw checkpoint after each completed cell.
    write_table(do.call(rbind, sim_list[seq_len(i)]), "task1_simulation_results")
    message(
      "  done T=",
      cell$T,
      " ",
      cell$method,
      " ",
      cell$boot,
      " b=",
      cell$block,
      ": P(qhat=5) = ",
      round(sum(qs == 5, na.rm = TRUE) / length(qs), 3),
      " (",
      sum(is.na(qs)),
      " failed)"
    )
  }
  sim_results <- do.call(rbind, sim_list)
  write_table(sim_results, "task1_simulation_results")
  sim_summary <- summarise_simulation(sim_results)
  rownames(sim_summary) <- NULL
  write_table(sim_summary, "task1_simulation_summary")
  # Paired SOBI vs AMUSE comparison on the shared datasets.
  paired_summary <- paired_method_comparison(sim_results)
  rownames(paired_summary) <- NULL
  write_table(paired_summary, "task1_paired_method_comparisons")
  
  ## ---- Part E: course data ---------------------------------------------------
  message("Part E: loading ", file.path(data_dir, "timeseries.dat"))
  Xts <- read_timeseries(file.path(data_dir, "timeseries.dat"))
  T_ts <- nrow(Xts)
  p_ts <- ncol(Xts)
  if (max(blocks) > T_ts)
    stop("Course-data block length exceeds the number of rows")
  write_table(whiten_check(Xts), "task1_application_whitening_check")
  descr <- data.frame(
    series = colnames(Xts),
    mean = colMeans(Xts),
    sd = apply(Xts, 2L, stats::sd)
  )
  write_table(descr, "task1_data_description")
  # Sensitivity grid for the final analysis.
  app_grid <- expand.grid(
    method = c("SOBI", "AMUSE"),
    boot = c("F", "S"),
    block = blocks,
    stringsAsFactors = FALSE
  )
  app_list <- vector("list", nrow(app_grid))
  for (i in seq_len(nrow(app_grid))) {
    cell <- app_grid[i, ]
    fit <- if (cell$boot == "F")
      ladle_estimate(
        Xts,
        as.character(cell$method),
        "F",
        block = cell$block,
        nboot = nboot_app,
        seed = 20261019 + i
      )
    else
      ladle_estimate(
        Xts,
        as.character(cell$method),
        "S",
        p_boot = 1 / cell$block,
        nboot = nboot_app,
        seed = 20261019 + i
      )
    app_list[i] <- list(fit)
    message(
      "  ",
      cell$method,
      " ",
      cell$boot,
      " b=",
      cell$block,
      ": qhat = ",
      fit$qhat,
      " (",
      fit$status,
      "; ",
      fit$n_boot_used,
      "/",
      fit$nboot,
      " bootstrap fits used; gap = ",
      signif(fit$g_gap, 3),
      ")"
    )
  }
  app_results <- do.call(rbind, lapply(seq_along(app_list), function(i) {
    fit <- app_list[[i]]
    data.frame(
      method = fit$method,
      boot = fit$boot,
      block = if (fit$boot == "F")
        fit$block
      else
        1 / fit$p_boot,
      qhat = fit$qhat,
      status = fit$status,
      nboot = fit$nboot,
      n_boot_used = fit$n_boot_used,
      n_boot_attempted = fit$n_boot_attempted,
      n_boot_failed = fit$n_boot_failed,
      seed = fit$seed,
      g_min = if (is.finite(fit$qhat))
        min(fit$g)
      else
        NA_real_,
      g_gap = fit$g_gap
    )
  }))
  write_table(app_results, "task1_application_results")
  app_criteria <- do.call(rbind, lapply(seq_along(app_list), function(i) {
    fit <- app_list[[i]]
    data.frame(
      method = fit$method,
      boot = fit$boot,
      block = app_grid$block[i],
      k = fit$k,
      f = fit$f,
      phi = fit$phi,
      g = fit$g,
      lambda2 = fit$lambda2,
      qhat = fit$qhat,
      status = fit$status
    )
  }))
  write_table(app_criteria, "task1_application_criteria")
  
  ## ---- save the completed analysis -------------------------------------------
  save(
    ladle_example,
    whiten_num,
    boot_demo,
    sim_results,
    sim_summary,
    paired_summary,
    app_results,
    app_list,
    app_grid,
    app_criteria,
    run_settings,
    Xdemo,
    Xts,
    descr,
    file = file.path(result_dir, "task1_sim.RData")
  )
  message("Task 1 completed. Project root: ", normalizePath(root))
  print(app_results, row.names = FALSE)
  invisible(app_list)
}

## ---- 7. supplementary check (no Monte Carlo rerun needed) -------------------
# Loads the completed results, adds a numerical check of the expected
# stationary-bootstrap block length and re-saves the results.
# The dimension boxplot is drawn separately in Section 10e below.
run_task1_extra <- function(root = locate_root()) {
  if (is.null(root) ||
      !dir.exists(root))
    stop("Could not find the project root")
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(root)
  for (d in c(table_dir, figure_dir, result_dir))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  env <- new.env()
  if (!file.exists(file.path(result_dir, "task1_sim.RData")))
    stop("Run run_task1() first: results/task1_sim.RData does not exist")
  load(file.path(result_dir, "task1_sim.RData"), envir = env)
  # Expected block length of the stationary bootstrap. Block lengths are
  # drawn from a geometric distribution with success probability p, so their
  # expectation is 1/p. Check numerically with 50,000 draws from Geom(1/20):
  set.seed(20261020)
  block_len_check <- data.frame(
    p = 1 / 20,
    theoretical = 20,
    empirical = mean(stats::rgeom(50000, prob = 1 / 20) + 1L)
  )
  utils::write.csv(
    block_len_check,
    file.path(table_dir, "task1_block_length_check.csv"),
    row.names = FALSE
  )
  assign("block_len_check", block_len_check, envir = env)
  save(
    list = ls(envir = env),
    file = file.path(result_dir, "task1_sim.RData"),
    envir = env
  )
  message("Task 1 extra checks done. block_len_check:")
  print(block_len_check, row.names = FALSE)
}

## ---- 8. entry point --------------------------------------------------------
# TRUE: load saved results and draw figures; FALSE: run the analysis first.
# Edit a figure and its own save_sci() call in Section 10, then rerun that
# figure's block. The data loaded in Section 9 remain available in the session.
plots_only <- FALSE
reduced_mc <- FALSE  # TRUE: reduced design; FALSE: original design
# Optional extra check: run run_task1_extra() manually when needed.
if (!plots_only) {
  # Set all four values so changing modes replaces previous MC settings.
  Sys.setenv(
    TASK1_NREP   = if (reduced_mc)
      "40"
    else
      "200",
    TASK1_NBOOT  = "100",
    TASK1_TS     = if (reduced_mc)
      "200,2000"
    else
      "200,500,1000,2000",
    TASK1_BLOCKS = if (reduced_mc)
      "20,40"
    else
      "10,20,40"
  )
  run_task1()
}

## ---- 9. load saved results for plotting ------------------------------------
# No plotting functions are defined: each figure below is ordinary code.
# Run this section once per session before editing individual figures.
# Captions use the settings saved with the completed analysis.
saved_results <- file.path(result_dir, "task1_sim.RData")
if (!file.exists(saved_results))
  stop("Saved results not found: set plots_only <- FALSE and run task1.R once")
task1_results <- new.env(parent = emptyenv())
load(saved_results, envir = task1_results)

# Add the supplementary result without rerunning the Monte Carlo study.
if (!exists("block_len_check", envir = task1_results, inherits = FALSE)) {
  run_task1_extra()
  load(saved_results, envir = task1_results)
}

boot_demo <- task1_results$boot_demo
ladle_example <- task1_results$ladle_example
sim_summary <- task1_results$sim_summary
sim_results <- task1_results$sim_results
app_criteria <- task1_results$app_criteria
run_settings <- task1_results$run_settings
nboot_demo <- run_settings$nboot_demo
nboot_example <- run_settings$nboot_example
nrep <- run_settings$nrep
nboot_screen <- run_settings$nboot_screen
T_ts <- nrow(task1_results$Xts)
p_ts <- ncol(task1_results$Xts)
q_true <- 5

## ---- 10a. Part B: bootstrap ACF comparison ----
resample_label <- if (is.null(nboot_demo)) {
  "Mean bootstrap ACF per scheme"
} else {
  paste0("Mean ACF over ", nboot_demo, " resamples per scheme")
}
pd <- do.call(rbind, lapply(c("original", "iid", "fixed_block", "stationary"), function(nm)
  data.frame(
    lag = boot_demo$lag,
    block = boot_demo$block,
    method = nm,
    acf = boot_demo[[nm]]
  )))
pd$method <- factor(pd$method,
                    levels = c("original", "iid", "fixed_block", "stationary"))
p_boot <- ggplot(pd, aes(lag, acf, colour = method, linetype = method)) +
  geom_hline(yintercept = 0,
                      colour = "grey70",
                      linewidth = .25) +
  geom_line(linewidth = .5) +
  facet_wrap(~ block, nrow = 1L, labeller = label_both) +
  labs(
    x = "Lag",
    y = "Autocorrelation",
    caption = paste0(
      resample_label,
      "; one Setting-1 signal, T = 500.\n",
      "For stationary bootstrap, block denotes its expected length."
    )
  ) +
  scale_colour_manual(
    values = c(
      original = "black",
      iid = PALETTE[2],
      fixed_block = PALETTE[1],
      stationary = PALETTE[3]
    )
  )

# Adjust this figure's width and height here.
save_sci(p_boot,
         "task1_fig_bootstrap_demo",
         width = "double",
         height_mm = 75)

## ---- 10b. Part C: ladle criterion example ----
p_ladle <- NULL
if (is.finite(ladle_example$qhat)) {
  lc <- do.call(rbind, lapply(c("f", "phi", "g"), function(nm)
    data.frame(
      k = ladle_example$k,
      criterion = nm,
      value = ladle_example[[nm]]
    )))
  p_ladle <- ggplot(lc,
                             aes(k, value, colour = criterion, shape = criterion)) +
    geom_vline(xintercept = q_true,
                        linetype = 2,
                        colour = "grey50") +
    geom_line(linewidth = .5) +
    geom_point(size = 1.6) +
    geom_point(
      data = data.frame(k = ladle_example$qhat, value = min(ladle_example$g)),
      aes(k, value),
      inherit.aes = FALSE,
      shape = 4,
      size = 3,
      stroke = 1
    ) +
    scale_x_continuous(breaks = ladle_example$k) +
    labs(
      x = "Component k",
      y = "Criterion",
      caption = paste0(
        "T = 1000, Setting 1; dashed line: true q = ",
        q_true,
        ".\n",
        "Cross: selected qhat = ",
        ladle_example$qhat,
        "; bootstrap fits used: ",
        ladle_example$n_boot_used,
        "/",
        nboot_example,
        "."
      )
    ) +
    scale_colour_manual(values = c(
      f = "black",
      phi = PALETTE[2],
      g = PALETTE[1]
    ))
  
  # Adjust this figure's width and height here.
  save_sci(p_ladle,
           "task1_fig_ladle_example",
           width = "double",
           height_mm = 75)
} else {
  message("Skipping the Part C figure: ", ladle_example$status)
}

## ---- 10c. Part D: Monte Carlo recovery ----
Ts <- sort(unique(sim_summary$T))
if (is.null(nrep))
  nrep <- sort(unique(sim_summary$n_attempted))
repetition_label <- if (length(nrep)) {
  paste0(paste(nrep, collapse = ", "), " Monte Carlo repetitions per cell")
} else {
  "Monte Carlo recovery over the saved configurations"
}
bootstrap_label <- if (is.null(nboot_screen)) {
  ""
} else {
  paste0("; ", nboot_screen, " bootstrap samples per fit")
}
sp <- sim_summary
sp$method_boot <- paste0(sp$method, " / ", sp$boot)
p_mc <- ggplot(sp, aes(factor(block), method_boot, fill = p_qhat5)) +
  geom_tile(colour = "white", linewidth = .4) +
  geom_text(aes(label = sprintf("%.2f", p_qhat5)), size = 3) +
  facet_wrap( ~ T, ncol = min(4L, length(Ts))) +
  scale_fill_gradient(low = "white",
                               high = PALETTE[1],
                               limits = c(0, 1)) +
  labs(
    x = "Block length",
    y = "Method / bootstrap",
    fill = "P(qhat = 5)",
    caption = paste0(
      repetition_label,
      bootstrap_label,
      ". True q = 5.\n",
      "Accuracy uses all attempts. MCSE, 95% Wilson intervals and failures are in the summary CSV."
    )
  ) +
  theme(legend.position = "bottom")

# Adjust this figure's width and height here.
save_sci(p_mc,
         "task1_fig_simulation",
         width = "double",
         height_mm = 75)

## ---- 10d. Part E: course-data criteria ----
data_label <- if (is.null(T_ts)) {
  paste0("Course data: p = ", p_ts)
} else {
  paste0("Course data: T = ", T_ts, ", p = ", p_ts)
}
ac <- do.call(rbind, lapply(c("f", "phi", "g"), function(nm) {
  data.frame(
    k = app_criteria$k,
    method = paste0(app_criteria$method, " / ", app_criteria$boot),
    block = app_criteria$block,
    criterion = nm,
    value = app_criteria[[nm]],
    qhat = app_criteria$qhat
  )
}))
ac$criterion <- factor(ac$criterion, levels = c("f", "phi", "g"))
ac_min <- ac[ac$criterion == "g" &
               is.finite(ac$qhat) &
               ac$k == ac$qhat, , drop = FALSE]
p_app <- ggplot(ac, aes(k, value, colour = method, linetype = method)) +
  geom_line(linewidth = .5) +
  geom_point(
    data = ac_min,
    shape = 4,
    size = 3,
    stroke = 1,
    na.rm = TRUE
  ) +
  facet_grid(criterion ~ block,
                      scales = "free_y",
                      labeller = label_both) +
  scale_x_continuous(breaks = seq(0, p_ts - 1L, by = 2L)) +
  labs(
    x = "Component k",
    y = "Criterion",
    caption = paste0(
      data_label,
      "; crosses mark the selected minima of g.\n",
      "F: fixed block; S: stationary bootstrap, with the displayed expected block length."
    )
  ) +
  scale_colour_manual(
    values = c(
      "SOBI / F" = PALETTE[1],
      "SOBI / S" = PALETTE[3],
      "AMUSE / F" = PALETTE[2],
      "AMUSE / S" = PALETTE[4]
    )
  )

# Adjust this figure's width and height here.
save_sci(p_app,
         "task1_fig_application",
         width = "double",
         height_mm = 120)

## ---- 10e. Part D: estimated-dimension distribution ----
sr <- sim_results
p_box <- ggplot(sr, aes(factor(T), qhat)) +
  geom_hline(yintercept = q_true,
                      linetype = 2,
                      colour = "grey50") +
  geom_boxplot(
    fill = PALETTE[1],
    alpha = .35,
    outlier.size = .4,
    linewidth = .3,
    na.rm = TRUE
  ) +
  facet_grid(method + boot ~ block, labeller = label_both) +
  labs(
    x = "Sample size T",
    y = "Estimated dimension qhat",
    caption = paste0(
      "Each panel has one method, bootstrap scheme and block length.\n",
      "Horizontal line: true q = ",
      q_true,
      "; failed fits excluded here and counted in the summary CSV."
    )
  )

# Adjust this figure's width and height here.
save_sci(p_box,
         "task1_fig_qhat_box",
         width = "double",
         height_mm = 120)

# To preview a figure in an interactive session, run its object name, e.g. p_mc.
# After changing its layers or labels, rerun that figure's code and save_sci().
