find_project_root <- function() {
  command_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  source_files <- unlist(lapply(sys.frames(), function(x) x$ofile))
  scripts <- c(source_files, sub("^--file=", "", command_file))
  scripts <- scripts[file.exists(scripts)]
  candidates <- c(
    if (length(scripts)) dirname(dirname(normalizePath(scripts))),
    getwd(), dirname(getwd())
  )
  valid <- vapply(candidates, function(path) {
    file.exists(file.path(path, "R", "task1_functions.R")) &&
      file.exists(file.path(path, "data", "timeseries.dat"))
  }, logical(1))
  if (!any(valid)) stop("Cannot locate the project directory.")
  normalizePath(candidates[which(valid)[1]])
}

project_root <- find_project_root()
source(file.path(project_root, "R", "task1_functions.R"))
rjd_iteration_limit <- 10000L
sos_decompose_part_c <- sos_decompose
sos_decompose <- function(C_list, rjd_eps = 1e-6,
                          rjd_maxiter = rjd_iteration_limit) {
  sos_decompose_part_c(C_list, rjd_eps = rjd_eps, rjd_maxiter = rjd_maxiter)
}
if (!requireNamespace("JADE", quietly = TRUE)) {
  stop("Install the JADE package before running this analysis.")
}
data_dir <- file.path(project_root, "output", "data")
table_dir <- file.path(project_root, "output", "tables")
figure_dir <- file.path(project_root, "output", "figures")
for (directory in c(data_dir, table_dir, figure_dir)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

read_integer_option <- function(name, default, minimum = 1L) {
  value <- suppressWarnings(as.numeric(Sys.getenv(name, as.character(default))))
  if (length(value) != 1L || !is.finite(value) ||
      value < minimum || value != floor(value) || value > .Machine$integer.max) {
    stop(name, " must be an integer >= ", minimum, ".")
  }
  as.integer(value)
}
n_repetitions <- read_integer_option("TASK1_REPS", 100L, 2L)
n_boot_simulation <- read_integer_option("TASK1_SIM_BOOT", 100L, 2L)
n_boot_data <- read_integer_option("TASK1_DATA_BOOT", 500L, 2L)
n_seed_runs <- read_integer_option("TASK1_SEED_RUNS", 5L, 2L)
n_cores <- read_integer_option("TASK1_CORES", 4L)
if (.Platform$OS.type == "windows") n_cores <- 1L
resume <- tolower(Sys.getenv("TASK1_RESUME", "true")) == "true"
sample_sizes <- c(200L, 500L, 1000L)
block_lengths <- c(10L, 20L, 40L)
methods <- c("fixed", "stationary")
lags <- 1:12
true_dimension <- 5L
primary_block_length <- 20L
master_seed <- 20261021L
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

write_task1_csv <- function(object, directory, filename) {
  stopifnot(startsWith(filename, "task1_"))
  write.csv(object, file.path(directory, filename), row.names = FALSE)
}

wilson_interval <- function(successes, trials, level = 0.95) {
  if (trials == 0L) return(c(lower = NA_real_, upper = NA_real_))
  z <- stats::qnorm(1 - (1 - level) / 2)
  probability <- successes / trials
  denominator <- 1 + z^2 / trials
  centre <- (probability + z^2 / (2 * trials)) / denominator
  half_width <- z * sqrt(probability * (1 - probability) / trials +
                           z^2 / (4 * trials^2)) / denominator
  c(lower = max(0, centre - half_width), upper = min(1, centre + half_width))
}

check_fit <- function(fit) {
  criterion <- fit$criterion
  p <- nrow(criterion)
  whitening_error <- max(abs(stats::cov(fit$whitening$whitened_data) - diag(p)))
  orthogonality_error <- max(abs(crossprod(fit$directions) - diag(p)))
  stopifnot(
    all(is.finite(as.matrix(criterion))),
    identical(criterion$dimension, 0:(p - 1L)),
    all(criterion$f >= -1e-12), all(criterion$phi >= -1e-12),
    max(abs(criterion$g - criterion$f - criterion$phi)) < 1e-12,
    fit$q_hat == criterion$dimension[which.min(criterion$g)],
    whitening_error < 1e-7, orthogonality_error < 1e-7
  )
  c(whitening_error = whitening_error, orthogonality_error = orthogonality_error)
}

run_checkpointed <- function(jobs, worker, config, filename, label) {
  path <- file.path(data_dir, filename)
  results <- vector("list", nrow(jobs))
  if (resume && file.exists(path)) {
    previous <- readRDS(path)
    if (is.null(previous$config$rjd_maxiter)) {
      previous$config$rjd_maxiter <- config$rjd_maxiter
    }
    if (identical(previous$config, config) && identical(previous$jobs, jobs)) {
      results <- previous$results
      failed_jobs <- vapply(results, function(result) {
        is.data.frame(result) && any(nzchar(result$error))
      }, logical(1))
      results[failed_jobs] <- list(NULL)
      message(label, ": resumed ", sum(!vapply(results, is.null, logical(1))), " jobs.")
    }
  }
  pending <- which(vapply(results, is.null, logical(1)))
  if (length(pending)) {
    batches <- split(pending, ceiling(seq_along(pending) / n_cores))
    for (batch in batches) {
      values <- parallel::mclapply(batch, function(i) worker(jobs[i, , drop = FALSE]),
                                  mc.cores = min(n_cores, length(batch)),
                                  mc.set.seed = FALSE)
      if (any(vapply(values, inherits, logical(1), "try-error"))) {
        failures <- values[vapply(values, inherits, logical(1), "try-error")]
        stop(label, ": a worker failed; completed batches remain in the checkpoint.\n",
             paste(unlist(failures), collapse = "\n"))
      }
      results[batch] <- values
      temporary <- paste0(path, ".tmp")
      saveRDS(list(config = config, jobs = jobs, results = results), temporary)
      if (!file.rename(temporary, path)) stop("Could not save checkpoint: ", path)
      message(label, ": ", sum(!vapply(results, is.null, logical(1))),
              "/", nrow(jobs), " jobs completed.")
    }
  }
  results
}

common_config <- list(
  analysis_version = 1L, lags = lags, block_lengths = block_lengths,
  methods = methods, master_seed = master_seed,
  functions_md5 = unname(tools::md5sum(file.path(project_root, "R", "task1_functions.R"))),
  R_version = R.version.string, JADE_version = as.character(utils::packageVersion("JADE"))
)

simulation_jobs <- expand.grid(
  replicate = seq_len(n_repetitions), sample_size = sample_sizes,
  KEEP.OUT.ATTRS = FALSE
)
simulation_config <- c(common_config, list(
  repetitions = n_repetitions, n_boot = n_boot_simulation,
  sample_sizes = sample_sizes, rjd_maxiter = rjd_iteration_limit
))

simulation_worker <- function(job) {
  size_id <- match(job$sample_size, sample_sizes)
  data_seed <- master_seed + size_id * 100000L + job$replicate
  X <- simulate_setting1(job$sample_size, seed = data_seed)
  rows <- vector("list", length(methods) * length(block_lengths))
  index <- 0L
  for (method in methods) {
    for (block_length in block_lengths) {
      index <- index + 1L
      bootstrap_seed <- master_seed + 1000000L + size_id * 1000000L +
        job$replicate * 10L + index
      started <- proc.time()[["elapsed"]]
      answer <- tryCatch({
        fit <- ladle_fit(X, lags, n_boot = n_boot_simulation,
                         bootstrap_method = method, block_length = block_length,
                         seed = bootstrap_seed)
        checks <- check_fit(fit)
        list(q_hat = fit$q_hat, error = "", checks = checks,
             criterion_gap = sort(fit$criterion$g)[2] - min(fit$criterion$g))
      }, error = function(e) {
        list(q_hat = NA_integer_, error = conditionMessage(e),
             checks = c(NA_real_, NA_real_), criterion_gap = NA_real_)
      })
      rows[[index]] <- data.frame(
        sample_size = job$sample_size, replicate = job$replicate,
        method = method, block_length = block_length, n_boot = n_boot_simulation,
        data_seed = data_seed, bootstrap_seed = bootstrap_seed,
        true_dimension = true_dimension, estimated_dimension = answer$q_hat,
        criterion_gap = answer$criterion_gap,
        whitening_error = answer$checks[1], orthogonality_error = answer$checks[2],
        elapsed_seconds = proc.time()[["elapsed"]] - started,
        error = answer$error, row.names = NULL
      )
    }
  }
  do.call(rbind, rows)
}

simulation_results <- do.call(rbind, run_checkpointed(
  simulation_jobs, simulation_worker, simulation_config,
  "task1_d_simulation_checkpoint.rds", "Part D"
))
write_task1_csv(simulation_results, data_dir, "task1_d_simulation_results.csv")
write_task1_csv(subset(simulation_results, nzchar(error)), table_dir,
                "task1_d_numerical_failures.csv")

simulation_groups <- split(simulation_results, interaction(
  simulation_results$sample_size, simulation_results$method,
  simulation_results$block_length, drop = TRUE
))
simulation_summary <- do.call(rbind, lapply(simulation_groups, function(group) {
  estimates <- group$estimated_dimension[!is.na(group$estimated_dimension)]
  n_valid <- length(estimates)
  n_correct <- sum(estimates == true_dimension)
  interval <- wilson_interval(n_correct, n_valid)
  probability <- if (n_valid) n_correct / n_valid else NA_real_
  data.frame(
    group[1, c("sample_size", "method", "block_length", "n_boot")],
    repetitions = nrow(group), successful_fits = n_valid,
    failed_fits = nrow(group) - n_valid,
    correct_rate = probability,
    correct_rate_all_attempts = n_correct / nrow(group),
    mcse = if (n_valid) sqrt(probability * (1 - probability) / n_valid) else NA_real_,
    ci_lower = interval[1], ci_upper = interval[2],
    under_rate = if (n_valid) mean(estimates < true_dimension) else NA_real_,
    over_rate = if (n_valid) mean(estimates > true_dimension) else NA_real_,
    mean_dimension = if (n_valid) mean(estimates) else NA_real_,
    sd_dimension = if (n_valid > 1) stats::sd(estimates) else NA_real_,
    mean_absolute_error = if (n_valid) mean(abs(estimates - true_dimension)) else NA_real_,
    row.names = NULL
  )
}))
rownames(simulation_summary) <- NULL
simulation_summary <- simulation_summary[order(simulation_summary$sample_size,
  simulation_summary$method, simulation_summary$block_length), ]
write_task1_csv(simulation_summary, table_dir, "task1_d_simulation_summary.csv")

dimension_frequencies <- do.call(rbind, lapply(simulation_groups, function(group) {
  counts <- tabulate(group$estimated_dimension + 1L, nbins = 10L)
  data.frame(
    sample_size = group$sample_size[1], method = group$method[1],
    block_length = group$block_length[1], estimated_dimension = 0:9,
    count = counts, proportion = if (sum(counts)) counts / sum(counts) else NA_real_
  )
}))
write_task1_csv(dimension_frequencies, table_dir, "task1_d_dimension_frequencies.csv")

paired_comparisons <- do.call(rbind, lapply(sample_sizes, function(n) {
  do.call(rbind, lapply(block_lengths, function(L) {
    group <- subset(simulation_results, sample_size == n & block_length == L)
    paired <- merge(subset(group, method == "fixed"),
                    subset(group, method == "stationary"), by = "replicate",
                    suffixes = c("_fixed", "_stationary"))
    paired <- paired[complete.cases(paired[, c("estimated_dimension_fixed",
                                               "estimated_dimension_stationary")]), ]
    difference <- as.integer(paired$estimated_dimension_stationary == true_dimension) -
      as.integer(paired$estimated_dimension_fixed == true_dimension)
    data.frame(sample_size = n, block_length = L, paired_repetitions = nrow(paired),
               stationary_minus_fixed_correct_rate = mean(difference),
               paired_mcse = stats::sd(difference) / sqrt(length(difference)),
               dimension_disagreement_rate = mean(paired$estimated_dimension_fixed !=
                                                   paired$estimated_dimension_stationary))
  }))
}))
write_task1_csv(paired_comparisons, table_dir, "task1_d_paired_method_comparison.csv")

method_colours <- c(fixed = "#2C7BB6", stationary = "#1A9641")
grDevices::png(file.path(figure_dir, "task1_d_recovery_probability.png"),
              width = 11, height = 4.3, units = "in", res = 220)
graphics::par(mfrow = c(1, 3), mar = c(4.3, 4.3, 3, 0.8))
for (L in block_lengths) {
  graphics::plot(NA, xlim = c(0.8, length(sample_sizes) + 0.2),
                 ylim = c(0, 1), xaxt = "n",
                 xlab = "Sample size n", ylab = "Pr(q-hat = 5)",
                 main = paste0("Block length = ", L))
  graphics::axis(1, at = seq_along(sample_sizes), labels = sample_sizes)
  graphics::abline(h = seq(0, 1, 0.2), col = "grey92")
  for (method in methods) {
    part <- simulation_summary[simulation_summary$block_length == L &
                                 simulation_summary$method == method, ]
    part <- part[order(part$sample_size), ]
    x <- seq_along(sample_sizes) + if (method == "fixed") -0.04 else 0.04
    graphics::arrows(x, part$ci_lower, x, part$ci_upper, angle = 90, code = 3,
                     length = 0.03, col = method_colours[method])
    graphics::lines(x, part$correct_rate, type = "o", pch = if (method == "fixed") 16 else 17,
                    col = method_colours[method])
  }
  graphics::legend("bottomright", c("Fixed", "Stationary", "95% Wilson intervals"),
                   col = c(method_colours, "grey45"), pch = c(16, 17, NA),
                   lty = c(1, 1, 1), bty = "n", cex = 0.75)
}
grDevices::dev.off()

dimension_colours <- c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6",
                       "#238B45", "#FDD0A2", "#FDAE6B", "#F16913", "#A63603")
grDevices::png(file.path(figure_dir, "task1_d_dimension_distribution.png"),
              width = 11, height = 4.4, units = "in", res = 220)
graphics::par(mfrow = c(1, 3), mar = c(4, 4, 3, 1), oma = c(2, 0, 0, 0))
for (n in sample_sizes) {
  values <- sapply(methods, function(method) {
    sapply(block_lengths, function(L) {
      part <- dimension_frequencies[dimension_frequencies$sample_size == n &
        dimension_frequencies$method == method & dimension_frequencies$block_length == L, ]
      part$proportion[order(part$estimated_dimension)]
    })
  }, simplify = "array")
  values <- matrix(values, nrow = 10)
  graphics::barplot(values, col = dimension_colours, border = NA,
                    names.arg = c(paste0("F", block_lengths), paste0("S", block_lengths)),
                    ylim = c(0, 1), ylab = "Proportion of repetitions", main = paste0("n = ", n))
}
graphics::par(fig = c(0, 1, 0, 0.06), new = TRUE, mar = rep(0, 4))
graphics::plot.new()
graphics::legend("center", legend = paste0("q=", 0:9), fill = dimension_colours,
                 horiz = TRUE, bty = "n", cex = 0.8)
grDevices::dev.off()

X <- as.matrix(read.table(file.path(project_root, "data", "timeseries.dat"), header = FALSE))
stopifnot(is.numeric(X), all(is.finite(X)), nrow(X) > max(lags), ncol(X) == 10L)
colnames(X) <- paste0("X", seq_len(ncol(X)))
application_jobs <- expand.grid(
  method = methods, block_length = block_lengths, seed_run = seq_len(n_seed_runs),
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
application_config <- c(common_config, list(
  n_boot = n_boot_data, seed_runs = n_seed_runs,
  data_md5 = unname(tools::md5sum(file.path(project_root, "data", "timeseries.dat"))),
  rjd_maxiter = rjd_iteration_limit
))
application_worker <- function(job) {
  seed <- master_seed + 10000000L + job$seed_run * 100L +
    match(job$method, methods) * 10L + match(job$block_length, block_lengths)
  fit <- ladle_fit(X, lags, n_boot = n_boot_data, bootstrap_method = job$method,
                   block_length = job$block_length, seed = seed)
  checks <- check_fit(fit)
  criterion <- fit$criterion
  ordering <- order(criterion$g)
  list(
    summary = data.frame(
      method = job$method, block_length = job$block_length, seed_run = job$seed_run,
      bootstrap_seed = seed, n = nrow(X), p = ncol(X), n_boot = n_boot_data,
      estimated_dimension = fit$q_hat, minimum_g = criterion$g[ordering[1]],
      runner_up_dimension = criterion$dimension[ordering[2]],
      criterion_gap = criterion$g[ordering[2]] - criterion$g[ordering[1]],
      whitening_error = checks[1], orthogonality_error = checks[2], row.names = NULL
    ),
    criterion = cbind(job, criterion),
    directions = if (job$seed_run == 1L) fit$directions else NULL,
    whitening = if (job$seed_run == 1L) fit$whitening else NULL,
    strengths = if (job$seed_run == 1L) fit$strengths else NULL
  )
}
application_results <- run_checkpointed(
  application_jobs, application_worker, application_config,
  "task1_e_application_checkpoint.rds", "Part E"
)
application_summary <- do.call(rbind, lapply(application_results, `[[`, "summary"))
application_criteria <- do.call(rbind, lapply(application_results, `[[`, "criterion"))
write_task1_csv(application_summary, table_dir, "task1_e_dimension_estimates.csv")
write_task1_csv(application_criteria, data_dir, "task1_e_ladle_criteria.csv")

application_stability <- do.call(rbind, lapply(split(application_summary, interaction(
  application_summary$method, application_summary$block_length, drop = TRUE
)), function(group) {
  counts <- table(group$estimated_dimension)
  data.frame(
    group[1, c("method", "block_length", "n_boot")], seed_runs = nrow(group),
    minimum_dimension = min(group$estimated_dimension),
    maximum_dimension = max(group$estimated_dimension),
    modal_dimension = as.integer(names(counts)[which.max(counts)]),
    modal_fraction = max(counts) / nrow(group),
    minimum_criterion_gap = min(group$criterion_gap),
    mean_criterion_gap = mean(group$criterion_gap), row.names = NULL
  )
}))
write_task1_csv(application_stability, table_dir, "task1_e_seed_stability.csv")

grDevices::png(file.path(figure_dir, "task1_e_ladle_components.png"),
              width = 11, height = 7, units = "in", res = 220)
graphics::par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
for (method in methods) {
  for (L in block_lengths) {
    part <- application_criteria[application_criteria$method == method &
      application_criteria$block_length == L & application_criteria$seed_run == 1L, ]
    q_hat <- part$dimension[which.min(part$g)]
    graphics::matplot(part$dimension, part[, c("f", "phi", "g")], type = "o",
                      pch = c(1, 2, 16), lty = c(2, 3, 1),
                      col = c("#2C7BB6", "#D7191C", "black"),
                      ylim = c(0, max(application_criteria$g) * 1.08),
                      xlab = "Candidate dimension k", ylab = "Criterion value",
                      main = paste0(method, ", L=", L, ": q-hat=", q_hat))
    graphics::abline(v = q_hat, col = "grey50", lty = 3)
    graphics::legend("topright", c("f(k)", "phi(k)", "g(k)"),
                     col = c("#2C7BB6", "#D7191C", "black"),
                     lty = c(2, 3, 1), pch = c(1, 2, 16), bty = "n", cex = 0.75)
  }
}
grDevices::dev.off()

grDevices::png(file.path(figure_dir, "task1_e_block_sensitivity.png"),
              width = 9, height = 4.6, units = "in", res = 220)
graphics::par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
length_colours <- c("#E69F00", "#0072B2", "#009E73")
for (method in methods) {
  graphics::plot(NA, xlim = c(0, 9), ylim = range(c(0, application_criteria$g)),
                 xlab = "Candidate dimension k", ylab = "g(k)", main = method)
  for (j in seq_along(block_lengths)) {
    part <- application_criteria[application_criteria$method == method &
      application_criteria$block_length == block_lengths[j], ]
    for (run in seq_len(n_seed_runs)) {
      one <- part[part$seed_run == run, ]
      graphics::lines(one$dimension, one$g,
                      col = grDevices::adjustcolor(length_colours[j], alpha.f = 0.25))
    }
    average <- aggregate(g ~ dimension, part, mean)
    graphics::lines(average$dimension, average$g, type = "o", pch = j,
                    lwd = 2, col = length_colours[j])
  }
  graphics::legend("topright", legend = paste0("L=", block_lengths),
                   col = length_colours, lty = 1, pch = 1:3, bty = "n", cex = 0.8)
}
grDevices::dev.off()

primary_id <- which(application_jobs$method == "fixed" &
  application_jobs$block_length == primary_block_length & application_jobs$seed_run == 1L)
primary <- application_results[[primary_id]]
sources <- primary$whitening$whitened_data %*% primary$directions
colnames(sources) <- paste0("component_", seq_len(ncol(sources)))
write_task1_csv(sources, data_dir, "task1_e_ordered_components.csv")
loadings <- primary$whitening$whitening_matrix %*% primary$directions
colnames(loadings) <- colnames(sources)
write_task1_csv(data.frame(variable = colnames(X), loadings), table_dir,
                "task1_e_unmixing_coefficients.csv")
component_summary <- data.frame(
  component = seq_len(ncol(X)), strength = primary$strengths,
  primary_classification = ifelse(seq_len(ncol(X)) <= primary$summary$estimated_dimension,
                                  "signal", "noise_candidate"),
  maximum_absolute_acf_lags_1_12 = vapply(seq_len(ncol(sources)), function(j) {
    max(abs(stats::acf(sources[, j], lag.max = 12, plot = FALSE)$acf[-1]))
  }, numeric(1)),
  lag_of_maximum_absolute_acf_1_12 = vapply(seq_len(ncol(sources)), function(j) {
    which.max(abs(stats::acf(sources[, j], lag.max = 12, plot = FALSE)$acf[-1]))
  }, integer(1))
)
write_task1_csv(component_summary, table_dir, "task1_e_component_summary.csv")
if (primary$summary$estimated_dimension > 0L) {
  last_signal <- component_summary[primary$summary$estimated_dimension, ]
  dominant_lag <- last_signal$lag_of_maximum_absolute_acf_1_12
  dependence_retention <- data.frame(
    component = last_signal$component, lag = dominant_lag,
    absolute_acf = last_signal$maximum_absolute_acf_lags_1_12,
    block_length = block_lengths,
    fixed_approximate_within_block_fraction = pmax(1 - dominant_lag / block_lengths, 0),
    stationary_no_restart_probability = (1 - 1 / block_lengths)^dominant_lag
  )
  write_task1_csv(dependence_retention, table_dir, "task1_e_dependence_retention.csv")
}
grDevices::png(file.path(figure_dir, "task1_e_component_acf.png"),
              width = 12, height = 6, units = "in", res = 220)
graphics::par(mfrow = c(2, 5), mar = c(4.2, 4.2, 3.2, 0.8),
              mgp = c(2.3, 0.7, 0))
for (j in seq_len(ncol(sources))) {
  stats::acf(sources[, j], lag.max = 24, ylim = c(-0.5, 1),
             main = paste0("Component ", j), ci = 0,
             col = if (j <= primary$summary$estimated_dimension) "#2C7BB6" else "grey40")
  graphics::abline(h = c(-1, 1) * 1.96 / sqrt(nrow(X)), col = "grey65", lty = 3)
}
grDevices::dev.off()

parameters <- data.frame(
  parameter = c("sample_sizes", "repetitions_per_sample_size", "simulation_bootstraps",
                "application_bootstraps", "application_seed_runs", "lags", "block_lengths",
                "primary_block_length", "master_seed", "cores", "rjd_maxiter", "R_version", "JADE_version"),
  value = c(paste(sample_sizes, collapse = ","), n_repetitions, n_boot_simulation,
            n_boot_data, n_seed_runs, paste(lags, collapse = ","),
            paste(block_lengths, collapse = ","), primary_block_length, master_seed,
            n_cores, rjd_iteration_limit, R.version.string,
            as.character(utils::packageVersion("JADE")))
)
write_task1_csv(parameters, table_dir, "task1_de_run_parameters.csv")
capture.output(utils::sessionInfo(), file = file.path(data_dir, "task1_de_session_info.txt"))
