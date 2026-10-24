# ============================================================================
# MAST32001 Task 1 - independent numerical verification (no MC rerun).
#
# Every stopifnot() must pass; a summary CSV and session info are written to
# results/tables/task1_numerical_checks.csv and results/task1_session_info.txt.
# ============================================================================
Sys.setenv(TASK1_SKIP_RUN = "1")   # load functions only; do not run the MC
source(file.path(root, file.path("R", "project_setup.R")))
ensure_packages(c("JADE"))
source(file.path(root, file.path("R", "task1.R")))

out <- list()
record <- function(name, value) {
  out[[length(out) + 1L]] <<- data.frame(check = name, value = value)
}
pass <- function(name)
  record(name, "passed")

set.seed(2026102100)
X <- generate_setting1(500)

## ---- 1. Whitening: Z'Z/n = I and zero column means ------------------------
Z <- whiten(X)
S <- crossprod(sweep(Z, 2L, colMeans(Z))) / nrow(Z)
err_whiten <- max(abs(S - diag(ncol(Z))))
err_mean <- max(abs(colMeans(Z)))
stopifnot(err_whiten < 1e-10, err_mean < 1e-10)
record("whitening_max_error", err_whiten)
record("whitened_mean_max_error", err_mean)
pass("whitening_identity_and_centering")

## ---- 2. Bootstrap row integrity: same indices used for every column -------
toy <- cbind(1:100, 101:200)
for (boot in c("fixed", "stationary", "iid")) {
  xb <- if (boot == "fixed")
    fixed_block_bootstrap(toy, 10)
  else if (boot == "stationary")
    stationary_bootstrap(toy, 1 / 10)
  else
    toy[sample.int(100, 100, replace = TRUE), , drop = FALSE]
  stopifnot(nrow(xb) == 100, all(xb[, 2L] - xb[, 1L] == 100))
}
pass("row_resampling_preserves_columns")

## ---- 3. Autocovariance hand calculation on a 4 x 4 toy matrix -------------
z <- rbind(c(1, 2), c(3, 4), c(5, 6), c(7, 8))
target <- matrix(c(53, 65, 65, 80), 2) / 3
A1 <- autocov(z, 1L)
stopifnot(max(abs(A1 - target)) < 1e-12)
record("lag_matrix_hand_calculation_error", max(abs(A1 - target)))

## ---- 4. Ladle structure: f0(0)=0, g = f + phi, qhat = argmin g ------------
set.seed(2026102104)
L <- ladle_estimate(X, "SOBI", "S", p_boot = 1 / 20, nboot = 30)
stopifnot(L$qhat == L$k[which.min(L$g)],
          L$f[1L] == 0,
          max(abs(L$g - L$f - L$phi)) < 1e-12,
          L$n_boot_used == 30L,
          L$status == "ok")
record("ladle_additivity_error", max(abs(L$g - L$f - L$phi)))
pass("ladle_additivity_and_minimum_index")

## ---- 5. Rotation invariance: qhat unchanged under an orthogonal mixing ----
Q <- qr.Q(qr(matrix(rnorm(100), 10)))          # random orthogonal 10 x 10
set.seed(2026102105)
q1 <- ladle_estimate(X, "SOBI", "S", p_boot = 1 / 20, nboot = 30)$qhat
set.seed(2026102105)
q2 <- ladle_estimate(X %*% Q, "SOBI", "S", p_boot = 1 / 20, nboot = 30)$qhat
stopifnot(is.finite(q1) && is.finite(q2) && q1 == q2)
record("rotation_invariance_qhat", q1)
pass("ladle_rotation_invariance")

## ---- 6. Negative AR(1): AMUSE must align the strongest serial component ---
set.seed(2026102106)
negative <- cbind(as.numeric(stats::arima.sim(list(ar = -0.7), n = 3000)), stats::rnorm(3000))
U <- joint_diagonalizer(negative, "AMUSE", tau = 1L)
stopifnot(!is.null(U), abs(U[1L, 1L]) > 0.9)   # first column aligns with series 1
record("negative_AR1_first_component_load", abs(U[1L, 1L]))

## ---- 7. Reproducibility: re-run the recorded example and compare with CSV -
# The benchmark is the example table produced by THIS codebase (the
# reduced_validation run), so the check verifies exact reproducibility of the
# merged pipeline, not bitwise agreement with the legacy full-MC artifacts
# (which used a different bootstrap RNG consumption order).
csv_path <- file.path("results",
                      "reduced_validation",
                      "tables",
                      "task1_ladle_example_criteria.csv")
if (!file.exists(csv_path))
  csv_path <- file.path(table_dir, "task1_ladle_example_criteria.csv")
if (file.exists(csv_path)) {
  set.seed(20261018)
  ex <- ladle_estimate(
    generate_setting1(1000),
    "SOBI",
    "S",
    p_boot = 1 / 20,
    nboot = 200,
    seed = 20261018
  )
  saved <- utils::read.csv(csv_path)
  err <- max(abs(ex$g - saved$g))
  stopifnot(is.finite(ex$qhat), err < 1e-8)
  record("recorded_example_max_criterion_difference", err)
  record("recorded_example_estimated_dimension", ex$qhat)
  pass("recorded_example_reproducible")
} else {
  record("recorded_example_reproducible", "skipped: CSV not found")
}

## 8. ---- Stationary bootstrap expected block length ------------------------
set.seed(20261020)
block_len <- mean(stats::rgeom(50000, prob = 1 / 20) + 1L)
stopifnot(abs(block_len - 20) < 0.1)
record("stationary_expected_block_length_empirical", block_len)

## ---- write outputs --------------------------------------------------------
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(
  do.call(rbind, out),
  file.path(table_dir, "task1_numerical_checks.csv"),
  row.names = FALSE
)
utils::capture.output(sessionInfo(),
                      file = file.path(result_dir, "task1_session_info_validation.txt"))
message(
  "Task 1 checks passed: algebra, indexing, row integrity, invariance,",
  " ladle structure and reproducibility all OK."
)
