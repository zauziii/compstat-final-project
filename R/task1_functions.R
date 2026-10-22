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

