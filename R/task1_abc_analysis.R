find_project_root <- function() {
  command_file <- grep(
    "^--file=",
    commandArgs(trailingOnly = FALSE),
    value = TRUE
  )

  candidates <- character()

  if (length(command_file) > 0) {
    script_path <- normalizePath(
      sub("^--file=", "", command_file[1]),
      mustWork = TRUE
    )
    candidates <- c(candidates, dirname(dirname(script_path)))
  }

  candidates <- c(candidates, getwd())

  candidates <- unique(normalizePath(candidates, mustWork = TRUE))
  valid <- vapply(
    candidates,
    function(path) {
      file.exists(file.path(path, "R", "task1_functions.R")) &&
        file.exists(file.path(path, "data", "timeseries.dat"))
    },
    logical(1)
  )

  if (!any(valid)) {
    stop("Run this script from the project root directory.")
  }

  candidates[which(valid)[1]]
}

project_root <- find_project_root()
source(file.path(project_root, "R","project_setup.R"))
source(file.path(project_root, "R", "task1_functions.R"))

data_dir <- file.path(project_root, "results", "data")
table_dir <- file.path(project_root, "results", "tables")
figure_dir <- file.path(project_root, "results", "figures")

dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

write_task1_csv <- function(object, directory, filename) {
  if (!startsWith(filename, "task1_")) {
    stop("Every output filename must begin with 'task1_'.")
  }

  write.csv(
    object,
    file = file.path(directory, filename),
    row.names = FALSE
  )
}

off_diagonal <- function(M) {
  M[row(M) != col(M)]
}

acf_vector <- function(x, lag_max = 20) {
  as.numeric(
    stats::acf(
      x,
      lag.max = lag_max,
      plot = FALSE,
      demean = TRUE
    )$acf
  )
}

mean_acf_error <- function(original, resampled, lag_max = 20) {
  original <- as.matrix(original)
  resampled <- as.matrix(resampled)

  variable_errors <- vapply(
    seq_len(ncol(original)),
    function(j) {
      original_acf <- acf_vector(original[, j], lag_max)[-1]
      resampled_acf <- acf_vector(resampled[, j], lag_max)[-1]
      mean(abs(original_acf - resampled_acf))
    },
    numeric(1)
  )

  mean(variable_errors)
}

joint_off_diagonal_ratio <- function(C_list, directions) {
  numerator <- sum(
    vapply(
      C_list,
      function(C) {
        rotated <- t(directions) %*% C %*% directions
        sum(off_diagonal(rotated)^2)
      },
      numeric(1)
    )
  )

  denominator <- sum(
    vapply(C_list, function(C) sum(C^2), numeric(1))
  )

  numerator / denominator
}

# Read the course data as 1,000 observations of 10 variables.
# The supplied file has no header row.
X <- as.matrix(
  read.table(
    file.path(project_root, "data", "timeseries.dat"),
    header = FALSE
  )
)
colnames(X) <- paste0("X", seq_len(ncol(X)))

# Part A: whitening and lagged autocovariance matrices
# ----------------------------------------------------
whitening <- whiten_ts(X)
Y <- whitening$whitened_data
colnames(Y) <- paste0("Y", seq_len(ncol(Y)))

S_Y <- stats::cov(Y)
identity_error <- S_Y - diag(ncol(Y))
lags <- 1:12
C_list <- autocov_set(Y, lags = lags, symmetrize = TRUE)

whitening_summary <- data.frame(
  n = nrow(X),
  p = ncol(X),
  maximum_absolute_whitened_mean = max(abs(colMeans(Y))),
  minimum_whitened_variance = min(diag(S_Y)),
  maximum_whitened_variance = max(diag(S_Y)),
  maximum_absolute_whitened_off_diagonal = max(abs(off_diagonal(S_Y))),
  maximum_absolute_identity_error = max(abs(identity_error))
)

autocovariance_summary <- data.frame(
  lag = lags,
  frobenius_norm = vapply(
    C_list,
    function(C) sqrt(sum(C^2)),
    numeric(1)
  ),
  maximum_absolute_entry = vapply(
    C_list,
    function(C) max(abs(C)),
    numeric(1)
  )
)

write_task1_csv(Y, data_dir, "task1_whitened_data.csv")
write_task1_csv(
  whitening_summary,
  table_dir,
  "task1_whitening_summary.csv"
)
write_task1_csv(
  autocovariance_summary,
  table_dir,
  "task1_autocovariance_summary.csv"
)

# Use the shared publication theme, palette and 600 dpi PNG exporter.
matrix_data <- function(M, panel) {
  grid <- expand.grid(row = seq_len(nrow(M)), column = seq_len(ncol(M)))
  transform(grid, value = as.vector(M), panel = panel)
}
covariance_plot_data <- rbind(
  matrix_data(stats::cor(X), "Original sample correlation"),
  matrix_data(S_Y, "Whitened sample covariance")
)
covariance_plot <- ggplot(covariance_plot_data, aes(row, column, fill = value)) +
  geom_tile() +
  facet_wrap(~panel, nrow = 1) +
  scale_fill_gradient2(low = plot_palette[1], mid = "white",
                       high = plot_palette[2], midpoint = 0, limits = c(-1, 1),
                       oob = scales::squish, name = "Value") +
  scale_x_continuous(breaks = seq_len(ncol(X)), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq_len(ncol(X)), expand = c(0, 0)) +
  coord_equal() + labs(x = "Variable", y = "Variable")
save_sci(covariance_plot, "task1_whitening_covariance", "double", height_mm = 75)

# Part B: bootstrap methods and preservation of local dependence
# ----------------------------------------------------------------
set.seed(20261013)
representative_variable <- which.max(
  vapply(
    seq_len(ncol(X)),
    function(j) abs(acf_vector(X[, j], lag_max = 1)[2]),
    numeric(1)
  )
)

iid_example <- iid_bootstrap(X)
fixed_example <- fixed_block_bootstrap(X, block_length = 20)
stationary_example <- stationary_bootstrap(X, mean_block_length = 20)

colnames(iid_example) <- colnames(X)
colnames(fixed_example) <- colnames(X)
colnames(stationary_example) <- colnames(X)

write_task1_csv(iid_example, data_dir, "task1_iid_bootstrap_example.csv")
write_task1_csv(
  fixed_example,
  data_dir,
  "task1_fixed_bootstrap_example.csv"
)
write_task1_csv(
  stationary_example,
  data_dir,
  "task1_stationary_bootstrap_example.csv"
)

block_lengths <- c(10, 20, 40)
n_dependence_boot <- 100
dependence_results <- vector(
  "list",
  1 + 2 * length(block_lengths)
)

dependence_results[[1]] <- data.frame(
  method = "iid",
  block_length = NA_integer_,
  replicate = seq_len(n_dependence_boot),
  mean_absolute_acf_error = vapply(
    seq_len(n_dependence_boot),
    function(b) mean_acf_error(X, iid_bootstrap(X)),
    numeric(1)
  )
)

position <- 2
for (block_length in block_lengths) {
  dependence_results[[position]] <- data.frame(
    method = "fixed",
    block_length = block_length,
    replicate = seq_len(n_dependence_boot),
    mean_absolute_acf_error = vapply(
      seq_len(n_dependence_boot),
      function(b) {
        X_boot <- fixed_block_bootstrap(
          X,
          block_length = block_length
        )
        mean_acf_error(X, X_boot)
      },
      numeric(1)
    )
  )
  position <- position + 1

  dependence_results[[position]] <- data.frame(
    method = "stationary",
    block_length = block_length,
    replicate = seq_len(n_dependence_boot),
    mean_absolute_acf_error = vapply(
      seq_len(n_dependence_boot),
      function(b) {
        X_boot <- stationary_bootstrap(
          X,
          mean_block_length = block_length
        )
        mean_acf_error(X, X_boot)
      },
      numeric(1)
    )
  )
  position <- position + 1
}

dependence_results <- do.call(rbind, dependence_results)
dependence_results$method_label <- ifelse(
  dependence_results$method == "iid",
  "iid",
  paste0(
    dependence_results$method,
    " (L=",
    dependence_results$block_length,
    ")"
  )
)

dependence_summary <- aggregate(
  mean_absolute_acf_error ~ method_label,
  data = dependence_results,
  FUN = function(x) {
    c(
      mean = mean(x),
      standard_deviation = stats::sd(x),
      median = stats::median(x)
    )
  }
)

dependence_summary <- data.frame(
  method = dependence_summary$method_label,
  mean_error = dependence_summary$mean_absolute_acf_error[, "mean"],
  standard_deviation = dependence_summary$mean_absolute_acf_error[
    , "standard_deviation"
  ],
  median_error = dependence_summary$mean_absolute_acf_error[, "median"],
  row.names = NULL
)

stationary_block_results <- do.call(
  rbind,
  lapply(
    block_lengths,
    function(block_length) {
      replicate_means <- vapply(
        seq_len(n_dependence_boot),
        function(b) {
          bootstrap_result <- stationary_bootstrap(
            X,
            mean_block_length = block_length,
            return_indices = TRUE
          )
          mean(bootstrap_result$observed_block_lengths)
        },
        numeric(1)
      )

      data.frame(
        target_expected_length = block_length,
        mean_observed_length = mean(replicate_means),
        standard_deviation = stats::sd(replicate_means)
      )
    }
  )
)

write_task1_csv(
  dependence_results,
  data_dir,
  "task1_bootstrap_dependence_results.csv"
)
write_task1_csv(
  dependence_summary,
  table_dir,
  "task1_bootstrap_dependence_summary.csv"
)
write_task1_csv(
  stationary_block_results,
  table_dir,
  "task1_stationary_block_summary.csv"
)

acf_lags <- 0:20
acf_comparison <- cbind(
  original = acf_vector(X[, representative_variable], 20),
  iid = acf_vector(iid_example[, representative_variable], 20),
  fixed = acf_vector(fixed_example[, representative_variable], 20),
  stationary = acf_vector(stationary_example[, representative_variable], 20)
)

acf_plot_data <- data.frame(lag = acf_lags, acf_comparison, check.names = FALSE) |>
  tidyr::pivot_longer(-lag, names_to = "method", values_to = "acf")
acf_plot_data$method <- factor(acf_plot_data$method,
                               levels = c("original", "iid", "fixed", "stationary"))
acf_labels <- c("Original", "iid", "Fixed L=20", "Stationary E[L]=20")
acf_plot <- ggplot(acf_plot_data, aes(lag, acf, colour = method,
                                      shape = method, linetype = method)) +
  geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.3) +
  geom_line(linewidth = 0.5) + geom_point(size = 1.6) +
  scale_colour_manual(values = c("#000000", plot_palette[c(2, 1, 3)]), labels = acf_labels) +
  scale_shape_manual(values = c(16, 1, 2, 5), labels = acf_labels) +
  scale_linetype_manual(values = c(1, 2, 3, 4), labels = acf_labels) +
  labs(x = "Lag", y = "Sample autocorrelation", colour = NULL, shape = NULL,
       linetype = NULL, title = paste0("Local dependence for variable X", representative_variable))
save_sci(acf_plot, "task1_bootstrap_acf_comparison", "double", height_mm = 80)

method_order <- c(
  "iid",
  "fixed (L=10)",
  "fixed (L=20)",
  "fixed (L=40)",
  "stationary (L=10)",
  "stationary (L=20)",
  "stationary (L=40)"
)
dependence_results$method_label <- factor(
  dependence_results$method_label,
  levels = method_order
)

error_plot <- ggplot(dependence_results,
                     aes(method_label, mean_absolute_acf_error, fill = method)) +
  geom_boxplot(linewidth = 0.3, outlier.size = 0.8) +
  scale_fill_manual(values = c(iid = plot_palette[2], fixed = plot_palette[1],
                               stationary = plot_palette[3]), guide = "none") +
  scale_x_discrete(labels = c("iid", "Fixed\nL=10", "Fixed\nL=20", "Fixed\nL=40",
                              "Stationary\nL=10", "Stationary\nL=20", "Stationary\nL=40")) +
  labs(x = NULL, y = "Mean absolute ACF error (lags 1–20)",
       title = "Preservation of local dependence")
save_sci(error_plot, "task1_bootstrap_acf_error", "double", height_mm = 75)

# Part C: f, phi and g checked on Setting 1 simulated data
# --------------------------------------------------------
X_simulated <- simulate_setting1(n = 1000, seed = 20261013)
write_task1_csv(
  X_simulated,
  data_dir,
  "task1_setting1_sample.csv"
)

fixed_ladle <- ladle_fit(
  X_simulated,
  lags = 1:12,
  n_boot = 100,
  bootstrap_method = "fixed",
  block_length = 20,
  seed = 101
)

stationary_ladle <- ladle_fit(
  X_simulated,
  lags = 1:12,
  n_boot = 100,
  bootstrap_method = "stationary",
  block_length = 20,
  seed = 102
)

fixed_criterion <- cbind(
  method = "fixed",
  fixed_ladle$criterion
)
stationary_criterion <- cbind(
  method = "stationary",
  stationary_ladle$criterion
)
ladle_criterion <- rbind(fixed_criterion, stationary_criterion)

part_c_checks <- data.frame(
  method = c("fixed", "stationary"),
  true_dimension = 5L,
  estimated_dimension = c(
    fixed_ladle$q_hat,
    stationary_ladle$q_hat
  ),
  whitening_identity_error = c(
    max(
      abs(
        stats::cov(fixed_ladle$whitening$whitened_data) -
          diag(ncol(X_simulated))
      )
    ),
    max(
      abs(
        stats::cov(stationary_ladle$whitening$whitened_data) -
          diag(ncol(X_simulated))
      )
    )
  ),
  direction_orthogonality_error = c(
    max(
      abs(
        crossprod(fixed_ladle$directions) - diag(ncol(X_simulated))
      )
    ),
    max(
      abs(
        crossprod(stationary_ladle$directions) - diag(ncol(X_simulated))
      )
    )
  ),
  joint_off_diagonal_ratio = c(
    joint_off_diagonal_ratio(
      fixed_ladle$autocovariances,
      fixed_ladle$directions
    ),
    joint_off_diagonal_ratio(
      stationary_ladle$autocovariances,
      stationary_ladle$directions
    )
  )
)

write_task1_csv(
  ladle_criterion,
  table_dir,
  "task1_ladle_criterion.csv"
)
write_task1_csv(
  part_c_checks,
  table_dir,
  "task1_part_c_checks.csv"
)

ladle_plot_data <- tidyr::pivot_longer(ladle_criterion, c(f, phi, g),
                                      names_to = "curve", values_to = "value")
ladle_panels <- setNames(paste0(c("Fixed", "Stationary"), ": q-hat = ",
                               part_c_checks$estimated_dimension), part_c_checks$method)
ladle_plot <- ggplot(ladle_plot_data, aes(dimension, value, colour = curve,
                                          shape = curve, linetype = curve)) +
  geom_vline(data = part_c_checks, aes(xintercept = estimated_dimension),
             colour = "grey50", linetype = 3, linewidth = 0.3) +
  geom_line(linewidth = 0.5) + geom_point(size = 1.7) +
  facet_wrap(~method, nrow = 1, labeller = as_labeller(ladle_panels)) +
  scale_colour_manual(values = c(f = plot_palette[1], phi = plot_palette[2], g = "#000000"),
                       breaks = c("f", "phi", "g"), labels = c("f(k)", "phi(k)", "g(k)")) +
  scale_shape_manual(values = c(f = 1, phi = 2, g = 16), breaks = c("f", "phi", "g"),
                      labels = c("f(k)", "phi(k)", "g(k)")) +
  scale_linetype_manual(values = c(f = 2, phi = 3, g = 1), breaks = c("f", "phi", "g"),
                         labels = c("f(k)", "phi(k)", "g(k)")) +
  scale_x_continuous(breaks = seq(0, 9, 2)) +
  labs(x = "Candidate dimension k", y = "Criterion value", colour = NULL,
       shape = NULL, linetype = NULL)
save_sci(ladle_plot, "task1_ladle_criterion", "double", height_mm = 75)

message("Task 1 A--C analysis completed.")
message("Data: ", data_dir)
message("Tables: ", table_dir)
message("Figures: ", figure_dir)
