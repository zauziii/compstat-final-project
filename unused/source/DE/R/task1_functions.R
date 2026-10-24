# Define a function for whitening a multivariate ts
#---------------------------------------------
whiten_ts <- function(X, tol = 1e-10) {
  
  X <- as.matrix(X)
  
  n <- nrow(X)
  p <- ncol(X)
  
  mu <- colMeans(X)
  X_centered <- sweep(X, 2, mu, FUN = "-")
  
  S <- crossprod(X_centered) / (n - 1)
  
  eig <- eigen(S, symmetric = TRUE)
  
  lambda <- eig$values
  V <- eig$vectors
  
  cutoff <- tol * max(lambda)
  
  if (any(lambda <= cutoff)) {
    stop(
      "The covariance matrix is singular or nearly singular."
    )
  }
  
  D_inv_sqrt <- diag(
    1 / sqrt(lambda),
    nrow = p,
    ncol = p
  )
  
  W <- V %*% D_inv_sqrt %*% t(V)
  
  Y <- X_centered %*% W
  
  list(
    mean = mu,
    covariance = S,
    whitening_matrix = W,
    whitened_data = Y
  )
}

# Define a function for estimating lagged autocovariance matrices
#---------------------------------------------
autocov_ts <- function(Y, lag, symmetrize = TRUE) {
  
  Y <- as.matrix(Y)
  n <- nrow(Y)
  
  if (lag < 1 || lag >= n) {
    stop("lag must satisfy 1 <= lag < n.")
  }
  
  Y_current <- Y[1:(n - lag), , drop = FALSE]
  Y_future <- Y[(lag + 1):n, , drop = FALSE]
  
  C_lag <- crossprod(Y_current, Y_future) / (n - lag)
  
  if (symmetrize) {
    C_lag <- (C_lag + t(C_lag)) / 2
  }
  
  C_lag
}

# Define a function for estimating autocovariance matrices at multiple lags
#---------------------------------------------
autocov_set <- function(Y, lags, symmetrize = TRUE) {
  
  Y <- as.matrix(Y)
  n <- nrow(Y)
  
  if (any(lags < 1) || any(lags >= n)) {
    stop("All lags must satisfy 1 <= lag < n.")
  }
  
  C_list <- lapply(
    lags,
    function(h) {
      autocov_ts(
        Y,
        lag = h,
        symmetrize = symmetrize
      )
    }
  )
  
  names(C_list) <- paste0("lag_", lags)
  
  return(C_list)
}

# Define an iid bootstrap function for comparison
# ------------------------------------------------
iid_bootstrap <- function(X, return_indices = FALSE) {
  
  X <- as.matrix(X)
  n <- nrow(X)

  indices <- sample.int(
    n,
    size = n,
    replace = TRUE
  )
  
  X_boot <- X[indices, , drop = FALSE]
  
  if (return_indices) {
    return(
      list(
        data = X_boot,
        indices = indices
      )
    )
  }
  
  return(X_boot)
}


# Define a fixed-block bootstrap function
# ------------------------------------------------
fixed_block_bootstrap <- function(
    X,
    block_length,
    return_indices = FALSE
) {
  
  X <- as.matrix(X)
  n <- nrow(X)
  
  if (
    length(block_length) != 1 ||
    block_length < 1 ||
    block_length > n ||
    block_length != floor(block_length)
  ) {
    stop(
      "block_length must be an integer between 1 and n."
    )
  }
  
  block_length <- as.integer(block_length)
  
  n_blocks <- ceiling(n / block_length)
  
  block_starts <- sample.int(
    n - block_length + 1,
    size = n_blocks,
    replace = TRUE
  )

  indices <- unlist(
    lapply(
      block_starts,
      function(start) {
        seq.int(
          from = start,
          length.out = block_length
        )
      }
    ),
    use.names = FALSE
  )
  
  indices <- indices[seq_len(n)]

  X_boot <- X[indices, , drop = FALSE]
  
  if (return_indices) {
    return(
      list(
        data = X_boot,
        indices = indices,
        block_starts = block_starts,
        block_length = block_length
      )
    )
  }
  
  return(X_boot)
}


# Define a stationary bootstrap function
# ------------------------------------------------
stationary_bootstrap <- function(
    X,
    mean_block_length,
    return_indices = FALSE
) {
  
  X <- as.matrix(X)
  n <- nrow(X)
  
  if (
    length(mean_block_length) != 1 ||
    mean_block_length < 1
  ) {
    stop("mean_block_length must be at least 1.")
  }
  
  restart_probability <- 1 / mean_block_length
  
  indices <- integer(n)
  
  new_block <- logical(n)
  
  indices[1] <- sample.int(n, size = 1)
  new_block[1] <- TRUE
  
  if (n >= 2) {
    for (t in 2:n) {
      
      if (runif(1) < restart_probability) {
        
        indices[t] <- sample.int(n, size = 1)
        new_block[t] <- TRUE
        
      } else {
        
        indices[t] <- if (
          indices[t - 1] == n
        ) {
          1L
        } else {
          indices[t - 1] + 1L
        }
      }
    }
  }
  
  X_boot <- X[indices, , drop = FALSE]
  
  if (return_indices) {
    
    block_start_positions <- which(new_block)
    
    block_end_positions <- c(
      block_start_positions[-1] - 1L,
      n
    )
    
    observed_block_lengths <-
      block_end_positions - block_start_positions + 1L
    
    return(
      list(
        data = X_boot,
        indices = indices,
        new_block = new_block,
        observed_block_lengths = observed_block_lengths,
        mean_block_length = mean_block_length,
        restart_probability = restart_probability
      )
    )
  }
  
  return(X_boot)
}


# Define a function for obtaining the SOS directions and component strengths
# ---------------------------------------------------------------------------
sos_decompose <- function(
    C_list,
    rjd_eps = 1e-6,
    rjd_maxiter = 10000
) {
  if (!is.list(C_list) || length(C_list) < 1) {
    stop("C_list must be a non-empty list of square matrices.")
  }

  C_list <- lapply(C_list, as.matrix)
  p <- nrow(C_list[[1]])

  valid_dimensions <- vapply(
    C_list,
    function(C) {
      nrow(C) == p && ncol(C) == p
    },
    logical(1)
  )

  if (!all(valid_dimensions)) {
    stop("All autocovariance matrices must have the same p by p dimensions.")
  }

  if (length(C_list) == 1) {
    diagonalization <- eigen(
      C_list[[1]],
      symmetric = TRUE
    )

    U <- diagonalization$vectors
  } else {
    if (!requireNamespace("JADE", quietly = TRUE)) {
      stop(
        "The JADE package is required for joint diagonalization."
      )
    }

    C_array <- simplify2array(C_list)

    diagonalization <- JADE::rjd(
      C_array,
      eps = rjd_eps,
      maxiter = rjd_maxiter
    )

    U <- diagonalization$V
  }

  rotated_diagonals <- vapply(
    C_list,
    function(C) {
      diag(t(U) %*% C %*% U)
    },
    numeric(p)
  )

  rotated_diagonals <- matrix(
    rotated_diagonals,
    nrow = p,
    ncol = length(C_list)
  )

  strengths <- rowSums(rotated_diagonals^2)
  ordering <- order(strengths, decreasing = TRUE)

  U <- U[, ordering, drop = FALSE]
  strengths <- strengths[ordering]
  rotated_diagonals <- rotated_diagonals[ordering, , drop = FALSE]

  list(
    directions = U,
    strengths = strengths,
    rotated_diagonals = rotated_diagonals,
    ordering = ordering
  )
}


# Define the eigenvalue-based part of the ladle criterion
# -------------------------------------------------------
compute_phi <- function(strengths) {
  strengths <- as.numeric(strengths)

  if (any(!is.finite(strengths)) || any(strengths < 0)) {
    stop("strengths must contain finite non-negative values.")
  }

  strengths <- sort(strengths, decreasing = TRUE)
  phi <- strengths / (1 + sum(strengths))

  return(phi)
}


# Define the bootstrap-variability part of the ladle criterion
# ------------------------------------------------------------
compute_f <- function(U_original, U_bootstrap) {
  U_original <- as.matrix(U_original)
  p <- ncol(U_original)

  if (!is.list(U_bootstrap) || length(U_bootstrap) < 1) {
    stop("U_bootstrap must be a non-empty list of direction matrices.")
  }

  valid_dimensions <- vapply(
    U_bootstrap,
    function(U) {
      U <- as.matrix(U)
      nrow(U) == nrow(U_original) && ncol(U) == p
    },
    logical(1)
  )

  if (!all(valid_dimensions)) {
    stop("All bootstrap direction matrices must match U_original.")
  }

  f0 <- numeric(p)
  f0[1] <- 0

  if (p >= 2) {
    for (k in seq_len(p - 1)) {
      B_original <- U_original[, seq_len(k), drop = FALSE]

      dissimilarities <- vapply(
        U_bootstrap,
        function(U_boot) {
          U_boot <- as.matrix(U_boot)
          B_boot <- U_boot[, seq_len(k), drop = FALSE]

          similarity <- abs(
            det(t(B_original) %*% B_boot)
          )

          similarity <- min(
            1,
            max(0, similarity)
          )

          1 - similarity
        },
        numeric(1)
      )

      f0[k + 1] <- mean(dissimilarities)
    }
  }

  f <- f0 / (1 + sum(f0))

  list(
    f0 = f0,
    f = f
  )
}


# Define the complete time-series ladle estimator
# ------------------------------------------------
ladle_fit <- function(
    X,
    lags,
    n_boot = 200,
    bootstrap_method = c("fixed", "stationary"),
    block_length = 20,
    symmetrize = TRUE,
    seed = NULL,
    progress = FALSE
) {
  X <- as.matrix(X)
  bootstrap_method <- match.arg(bootstrap_method)

  if (n_boot < 1 || n_boot != floor(n_boot)) {
    stop("n_boot must be a positive integer.")
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  original_whitening <- whiten_ts(X)
  Y <- original_whitening$whitened_data

  original_autocovariances <- autocov_set(
    Y,
    lags = lags,
    symmetrize = symmetrize
  )

  original_decomposition <- sos_decompose(
    original_autocovariances
  )

  U_original <- original_decomposition$directions
  phi <- compute_phi(
    original_decomposition$strengths
  )

  U_bootstrap <- vector(
    mode = "list",
    length = n_boot
  )

  for (b in seq_len(n_boot)) {
    X_boot <- switch(
      bootstrap_method,
      fixed = fixed_block_bootstrap(
        X,
        block_length = block_length
      ),
      stationary = stationary_bootstrap(
        X,
        mean_block_length = block_length
      )
    )

    bootstrap_whitening <- whiten_ts(X_boot)

    bootstrap_autocovariances <- autocov_set(
      bootstrap_whitening$whitened_data,
      lags = lags,
      symmetrize = symmetrize
    )

    bootstrap_decomposition <- sos_decompose(
      bootstrap_autocovariances
    )

    U_bootstrap[[b]] <- bootstrap_decomposition$directions

    if (
      progress &&
      (b %% max(1, floor(n_boot / 10)) == 0 || b == n_boot)
    ) {
      message(
        "Completed ",
        b,
        " of ",
        n_boot,
        " bootstrap samples."
      )
    }
  }

  f_result <- compute_f(
    U_original,
    U_bootstrap
  )

  f <- f_result$f
  g <- f + phi
  q_hat <- which.min(g) - 1L

  criterion <- data.frame(
    dimension = 0:(ncol(X) - 1),
    f = f,
    phi = phi,
    g = g
  )

  list(
    q_hat = q_hat,
    criterion = criterion,
    f0 = f_result$f0,
    directions = U_original,
    strengths = original_decomposition$strengths,
    whitening = original_whitening,
    autocovariances = original_autocovariances,
    bootstrap_directions = U_bootstrap,
    bootstrap_method = bootstrap_method,
    block_length = block_length,
    lags = lags,
    n_boot = n_boot
  )
}


# Generate one sample from Setting 1 of Nordhausen and Virta (2018)
# -----------------------------------------------------------------
simulate_setting1 <- function(n, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (n < 10 || n != floor(n)) {
    stop("n must be an integer of at least 10.")
  }

  signals <- cbind(
    as.numeric(
      arima.sim(
        model = list(
          ar = 0.6,
          ma = c(0, -0.4)
        ),
        n = n
      )
    ),
    as.numeric(
      arima.sim(
        model = list(
          ar = c(0.4, 0.1, 0.3),
          ma = c(0.2, 0.4)
        ),
        n = n
      )
    ),
    as.numeric(
      arima.sim(
        model = list(ar = 0.4),
        n = n
      )
    ),
    as.numeric(
      arima.sim(
        model = list(
          ma = c(0.5, 0.3, 0.1, 0.2)
        ),
        n = n
      )
    ),
    as.numeric(
      arima.sim(
        model = list(
          ar = c(0.1, 0.1, 0.3)
        ),
        n = n
      )
    )
  )

  noise <- matrix(
    rnorm(n * 5),
    nrow = n,
    ncol = 5
  )

  X <- cbind(signals, noise)
  X <- scale(X, center = TRUE, scale = TRUE)

  colnames(X) <- c(
    paste0("signal_", 1:5),
    paste0("noise_", 1:5)
  )

  return(as.matrix(X))
}

