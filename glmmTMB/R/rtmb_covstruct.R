## RTMB random-effect likelihoods, simulation, and covariance reporting.

## Partition the concatenated random effects and covariance parameters by term
## Term slicing is translated from allterms_nll() in glmmTMB.cpp:803-826
allterms_nll <- function(u, theta, terms) {
  "[<-" <- RTMB::ADoverload("[<-")

  nll <- 0
  corr <- vector("list", length(terms))
  sd <- vector("list", length(terms))
  fact_load <- vector("list", length(terms))
  names(corr) <- names(terms)
  names(sd) <- names(terms)
  names(fact_load) <- names(terms)

  if (length(terms) == 0) {
    output_u <- if (inherits(u, "simref")) u$value else u
    return(list(
      nll = nll, corr = corr, sd = sd, fact_load = fact_load, u = output_u
    ))
  }

  transformed_u <- u
  upointer <- 0L
  tpointer <- 0L
  np <- 0L

  for (i in seq_along(terms)) {
    term <- terms[[i]]
    nr <- term$blockSize * term$blockReps
    ## A zero-length theta block reuses the prev term's covariance parameters
    emptyTheta <- term$blockNumTheta == 0

    if (!emptyTheta) {
      np <- term$blockNumTheta
      theta_start <- tpointer + 1L
    } else {
      theta_start <- tpointer - np + 1L
    }

    useg <- u[(upointer + 1L):(upointer + nr)]

    if (np > 0) {
        tseg <- theta[theta_start:(theta_start + np - 1L)]
    } else {
      tseg <- numeric(0)
    }

    ans <- termwise_nll(useg, tseg, term)
    nll <- nll + ans$nll
    corr[[i]] <- ans$corr
    sd[[i]] <- ans$sd
    fact_load[[i]] <- ans$fact_load
    if (!inherits(transformed_u, "simref")) {
      transformed_u[(upointer + 1L):(upointer + nr)] <- ans$u
    }

    upointer <- upointer + nr
    tpointer <- tpointer + term$blockNumTheta
  }

  output_u <- if (inherits(transformed_u, "simref")) {
    transformed_u$value
  } else {
    transformed_u
  }

  list(
    nll = nll,
    corr = corr,
    sd = sd,
    fact_load = fact_load,
    u = output_u
  )
}

## Construct the correlation matrix used by TMB's
## density::UNSTRUCTURED_CORR_t. TMB fills the lower triangle row-wise,
## whereas matrix lower-triangle assignment in R fills it column-wise
## leading to a different theta ordering for dim >= 4
tmb_unstructured_corr <- function(n, theta) {
  expected <- n * (n - 1L) / 2L
  if (length(theta) != expected) {
    stop(
      "Expected ", expected, " correlation parameters for unstructured ",
      n, " by ", n, " correlation matrix, got ", length(theta)
    )
  }

  lower_idx <- which(lower.tri(diag(n)), arr.ind = TRUE)
  rowwise_idx <- lower_idx[order(lower_idx[, 1L], lower_idx[, 2L]),
                           , drop = FALSE]
  tmb_to_rtmb_order <- match(
    paste(lower_idx[, 1L], lower_idx[, 2L]),
    paste(rowwise_idx[, 1L], rowwise_idx[, 2L])
  )
  RTMB::unstructured(n)$corr(theta[tmb_to_rtmb_order])
}

## Simulation factor for density::UNSTRUCTURED_CORR_t, used in
## glmmTMB.cpp:407-425. This preserves the C++ simulation order for us()
## random effects instead of delegating simulation to RTMB::dmvnorm().
tmb_unstructured_sim_factor <- function(n, theta) {
  expected <- n * (n - 1L) / 2L
  if (length(theta) != expected) {
    stop(
      "Expected ", expected, " correlation parameters for unstructured ",
      n, " by ", n, " simulation factor, got ", length(theta)
    )
  }

  ## TMB fills the lower triangle row-wise.  Assigning theta to the upper
  ## triangle column-wise and transposing gives the same ordering, while direct
  ## lower.tri() assignment would use R's column-wise lower-triangle order.
  L <- diag(n)
  L[upper.tri(L)] <- theta
  L <- t(L)

  L / sqrt(rowSums(L * L))
}

simulate_tmb_unstructured <- function(U, sd, corr_par) {
  n <- nrow(U)
  reps <- ncol(U)
  sim_factor <- tmb_unstructured_sim_factor(n, as.vector(corr_par))
  sim_sd <- as.vector(sd)

  for (j in seq_len(reps)) {
    U_col <- U[, j]
    U_col[] <- sim_sd * as.vector(sim_factor %*% stats::rnorm(n))
  }
}

## Evaluate one random-effects term under its covariance structure
## Translation of the currently supported cases in
## termwise_nll(), glmmTMB.cpp:358-799
termwise_nll <- function(U, theta, term) {
  ## Preserve automatic differentiation when filling correlation matrices
  "[<-" <- RTMB::ADoverload("[<-")

  block_code <- term$blockCode
  name <- term$blockName %||% {
    block_name <- names(block_code)
    if (length(block_name) == 0L) {
      names(.valid_covstruct)[match(block_code, .valid_covstruct)]
    } else {
      block_name[1L]
    }
  }
  supported <- c(
    "diag", "homdiag", "us", "cs", "homcs", "toep", "homtoep",
    "ar1", "hetar1", "ou", "exp", "gau", "mat", "rr", "propto", "equalto"
  )

  if (!name %in% supported) {
    stop(
      "covariance structure not yet implemented: ", name,
      "; implemented structures are: ", paste(supported, collapse = ", ")
    )
  }

  n <- term$blockSize
  reps <- term$blockReps
  dim(U) <- c(n, reps)

  rr_rank <- NA_integer_
  if (name == "rr") {
    ntheta <- length(theta)
    rank_discriminant <- (2 * n + 1)^2 - 8 * ntheta
    if (rank_discriminant < 0) {
      stop(
        "Invalid covariance parameter count for 'rr': ", ntheta,
        "; rank discriminant is ", rank_discriminant,
        ", so no real-valued rank can be inferred for block size ", n
      )
    }
    rank_value <- (
      2 * n + 1 - sqrt(rank_discriminant)
    ) / 2
    rr_rank <- as.integer(round(rank_value))
    valid_rank <- is.finite(rank_value) &&
      abs(rank_value - rr_rank) < sqrt(.Machine$double.eps) &&
      rr_rank >= 1L &&
      rr_rank <= n

    if (!valid_rank) {
      stop(
        "Invalid covariance parameter count for 'rr': ", ntheta,
        "; inferred rank value is ", rank_value,
        ", rounded rank is ", rr_rank,
        ", valid ranks are integers from 1 to ", n
      )
    }
  }

  expected_num_theta <- switch(
    name,
    diag = n,
    homdiag = 1L,
    us = n * (n + 1L) / 2L,
    cs = n + 1L,
    homcs = 2L,
    toep = 2L * n - 1L,
    homtoep = n,
    ar1 = 2L,
    hetar1 = n + 1L,
    ou = 2L,
    exp = 2L,
    gau = 2L,
    mat = 3L,
    rr = n * rr_rank - (rr_rank - 1L) * rr_rank / 2L,
    propto = n * (n + 1L) / 2L + 1L,
    equalto = n * (n + 1L) / 2L
  )
  if (length(theta) != expected_num_theta) {
    stop(
      "Expected ", expected_num_theta, " covariance parameters for '",
      name, "', got ", length(theta)
    )
  }

  if (name == "rr") {
    ## Reduced-rank covariance; glmmTMB.cpp:698-761. The optimized random
    ## effects are spherical, while the linear predictor uses Lambda %*% u.
    nll <- 0
    simulation <- inherits(U, "simref")

    if (simulation && !term$simCode %in% .valid_simcode) {
      stop(
        "unknown simCode for rr covariance structure: ", term$simCode,
        "; known simCodes are: ",
        paste(names(.valid_simcode), .valid_simcode, sep = "=",
              collapse = ", ")
      )
    }

    if (!simulation || term$simCode == .valid_simcode[["random"]]) {
      for (j in seq_len(reps)) {
        nll <- nll - sum(RTMB::dnorm(U[, j], 0, 1, log = TRUE))
      }
    } else if (term$simCode == .valid_simcode[["zero"]]) {
      U[] <- 0
    } else {
      U[] <- U$getOrig(seq_along(U))
    }

    Lambda <- matrix(0, n, rr_rank)
    lam_diag <- head(theta, rr_rank)
    lam_lower <- utils::tail(theta, length(theta) - rr_rank)

    Lambda[row(Lambda) == col(Lambda)] <- lam_diag
    Lambda[row(Lambda) > col(Lambda)] <- lam_lower

    if (term$simCode != .valid_simcode[["fix"]]) {
      for (j in seq_len(reps)) {
        transformed_column <- Lambda %*% U[seq_len(rr_rank), j]
        if (simulation) {
          U_column <- U[, j]
          U_column[] <- transformed_column
        } else {
          U[, j] <- transformed_column
        }
      }
    }

    report_corr <- matrix(numeric(0), 0, 0)
    report_sd <- numeric(0)
    if (term$fullCor == 1L) {
      covariance <- Lambda %*% t(Lambda)
      report_sd <- sqrt(diag(covariance))
      report_corr <- covariance /
        (report_sd %*% t(report_sd))
    }

    return(list(
      nll = nll,
      corr = report_corr,
      sd = report_sd,
      fact_load = Lambda,
      u = if (simulation) NULL else as.vector(U)
    ))
  }

  ## Homogeneous structures use one standard-deviation parameter;
  ## heterogeneous structures use one parameter per term component.
  homogeneous <- c(
    "homdiag", "homcs", "homtoep", "ar1", "ou", "exp", "gau", "mat"
  )
  hetvar <- !name %in% homogeneous
  n_sd_par <- if (hetvar) n else 1L

  logsd <- if (hetvar) {
    head(theta, n)
  } else {
    rep(theta[1L], n)
  }

  sd <- exp(logsd)
  corr_par <- theta[-seq_len(n_sd_par)]

  ## propto uses an unstructured correlation matrix with an additional
  ## parameter that proportionally scales the covariance matrix.
  if (name == "propto") {
    loglambda <- utils::tail(corr_par, 1L)
    corr_par <- head(corr_par, -1L)
    sd <- exp(logsd + loglambda / 2)
  }

  ## Remove the "hom" prefix because homogeneous and heterogeneous
  ## variants differ only in their standard-deviation parameterization.
  cov_structure <- sub("^hom", "", name)

  ## propto and equalto use the unstructured correlation parameterization.
  density_structure <- if (cov_structure %in% c("propto", "equalto")) {
    "us"
  } else if (cov_structure == "hetar1") {
    "ar1"
  } else {
    cov_structure
  }

  C <- switch(
    density_structure,

    ## Diagonal covariance; glmmTMB.cpp:358-405
    diag = {
      ## Empty matrix means the diagonal structure has no correlation matrix;
      ## downstream code expects a matrix-valued placeholder rather than NULL.
      matrix(numeric(0), 0, 0)
    },

    ## Unstructured covariance; glmmTMB.cpp:407-440
    us = {
      tmb_unstructured_corr(n, corr_par)
    },

    ## Compound-symmetry covariance; glmmTMB.cpp:441-473
    cs = {
      a <- 1 / (n - 1)
      rho <- (1 / (1 + exp(-corr_par[1L]))) * (1 + a) - a
      corr <- diag(n)
      corr[row(corr) != col(corr)] <- rho
      corr
    },

    ## Toeplitz covariance; glmmTMB.cpp:474-506
    toep = {
      corr_params <- corr_par / sqrt(1 + corr_par^2)
      lag <- abs(row(diag(n)) - col(diag(n)))
      corr <- diag(n)
      off_diagonal <- lag > 0
      corr[off_diagonal] <- corr_params[lag[off_diagonal]]
      corr
    },

    ## Homogeneous AR(1) covariance; glmmTMB.cpp:507-590
    ar1 = {
      phi <- corr_par[1L] / sqrt(1 + corr_par[1L]^2)
      matrix(numeric(0), 0, 0)
    },

    ## OU covariance; glmmTMB.cpp:593-650
    ou = {
      times <- term$times
      if (length(times) != n) {
        stop(
          "OU time vector length must equal block size; got length(times)=",
          length(times), " and blockSize=", n
        )
      }
      decay <- exp(corr_par[1L])
      if (term$fullCor == 0) {
        matrix(numeric(0), 0, 0)
      } else {
        time_dist <- abs(outer(times, times, "-"))
        exp(-decay * time_dist)
      }
    },

    ## Exponential spatial covariance; glmmTMB.cpp:653-700
    exp = {
      spatial_dist <- term$dist
      spatial_dim <- dim(spatial_dist)
      if (length(spatial_dim) != 2L || any(spatial_dim != n)) {
        stop(
          "Dimension of distance matrix must equal block size for ", name,
          "; got dim(dist)=",
          paste(spatial_dim, collapse = " x "),
          " and blockSize=", n
        )
      }
      corr <- exp(-spatial_dist * exp(-corr_par[1L]))
      corr
    },

    ## Gaussian spatial covariance; glmmTMB.cpp:653-700
    gau = {
      spatial_dist <- term$dist
      spatial_dim <- dim(spatial_dist)
      if (length(spatial_dim) != 2L || any(spatial_dim != n)) {
        stop(
          "Dimension of distance matrix must equal block size for ", name,
          "; got dim(dist)=",
          paste(spatial_dim, collapse = " x "),
          " and blockSize=", n
        )
      }
      corr <- exp(-(spatial_dist^2) * exp(-2 * corr_par[1L]))
      corr
    },

    ## Matern covariance; glmmTMB.cpp:653-700
    mat = {
      spatial_dist <- term$dist
      spatial_dim <- dim(spatial_dist)
      if (length(spatial_dim) != 2L || any(spatial_dim != n)) {
        stop(
          "Dimension of distance matrix must equal block size for ", name,
          "; got dim(dist)=",
          paste(spatial_dim, collapse = " x "),
          " and blockSize=", n
        )
      }
      range <- exp(corr_par[1L])
      smoothness <- exp(corr_par[2L])
      scaled_dist <- spatial_dist / range
      diagonal <- row(scaled_dist) == col(scaled_dist)
      scaled_dist[diagonal] <- 1
      matern_corr <- RTMB::Vectorize(
        function(d) {
          d^smoothness * RTMB::besselK(d, smoothness) /
            (exp(lgamma(smoothness)) * 2^(smoothness - 1))
        },
        vectorize.args = "d"
      )
      corr <- matern_corr(as.vector(scaled_dist))
      dim(corr) <- c(n, n)
      corr[diagonal] <- 1
      corr
    },
    stop(
      "covariance structure not yet implemented: ", name,
      "; implemented density structures are: diag, us, cs, toep, ar1, ",
      "ou, exp, gau, mat"
    )
  )

  simulation <- inherits(U, "simref")
  simulate_density <- TRUE
  if (simulation) {
    if (!term$simCode %in% .valid_simcode) {
      stop(
        "unknown simCode for ", name, " covariance structure: ",
        term$simCode,
        "; known simCodes are: ",
        paste(names(.valid_simcode), .valid_simcode, sep = "=",
              collapse = ", ")
      )
    }

    ## Only the first group currently has explicit RTMB implementations for all
    ## C++ simulation modes.  The second group can simulate new random effects,
    ## but fixed/zero simulation would need separate structure-specific code.
    flexible_simulation <- c("diag", "us", "ar1", "hetar1", "ou")
    random_only_simulation <- c(
      "homdiag", "cs", "homcs", "toep", "homtoep",
      "exp", "gau", "mat"
    )

    if (name %in% flexible_simulation) {
      if (term$simCode == .valid_simcode[["zero"]]) {
        U[] <- 0
        simulate_density <- FALSE
      } else if (term$simCode == .valid_simcode[["fix"]]) {
        U[] <- U$getOrig(seq_along(U))
        simulate_density <- FALSE
      } else if (
        name == "us" && term$simCode == .valid_simcode[["random"]]
      ) {
        simulate_tmb_unstructured(U, sd, corr_par)
        simulate_density <- FALSE
      }
    } else if (
      name %in% random_only_simulation &&
      term$simCode != .valid_simcode[["random"]]
    ) {
      stop(
        "simCode '",
        names(.valid_simcode)[match(term$simCode, .valid_simcode)],
        "' is not implemented for ", name,
        " covariance structure; only random simulation is currently supported"
      )
    }
  }

  ## Diagonal structures factor into univariate normal densities;
  ## correlated structures use a scaled multivariate normal density.
  if (!simulate_density) {
    nll <- 0
  } else if (density_structure == "diag") {
    nll <- -sum(RTMB::dnorm(U, 0, sd, log = TRUE))
  } else if (density_structure == "ar1") {
    ## Match the state-space AR1 likelihood used by glmmTMB.cpp:522-554.
    ## This avoids building a dense covariance matrix and avoids routing the
    ## homogeneous AR1 case through RTMB::dautoreg()'s vector-scale wrapper.
    nll <- 0
    innovation_sd <- sqrt(1 - phi * phi)

    for (k in seq_len(reps)) {
      if (simulation) {
        U_sim <- numeric(n)
        if (name == "hetar1") {
          U_sim[1L] <- sd[1L] * stats::rnorm(1L)
          if (n > 1L) {
            for (j in 2:n) {
              U_sim[j] <- sd[j] * stats::rnorm(
                1L,
                mean = phi * U_sim[j - 1L] / sd[j - 1L],
                sd = innovation_sd
              )
            }
          }
        } else {
          U_sim[1L] <- stats::rnorm(1L, mean = 0, sd = sd[1L])
          if (n > 1L) {
            for (j in 2:n) {
              U_sim[j] <- stats::rnorm(
                1L,
                mean = phi * U_sim[j - 1L],
                sd = sd[1L] * innovation_sd
              )
            }
          }
        }
        U_column <- U[, k]
        U_column[] <- U_sim
      }
    }
    if (!simulation && name == "hetar1") {
      nll <- -sum(RTMB::dnorm(U[1L, ] / sd[1L], 0, 1, log = TRUE)) +
        reps * logsd[1L]
      if (n > 1L) {
        nll <- nll - sum(RTMB::dnorm(
          U[-1L, , drop = FALSE] / sd[-1L],
          phi * U[-n, , drop = FALSE] / sd[-n],
          innovation_sd,
          log = TRUE
        )) + reps * sum(logsd[-1L])
      }
    } else if (!simulation) {
      nll <- -sum(RTMB::dnorm(U[1L, ], 0, sd[1L], log = TRUE))
      if (n > 1L) {
        nll <- nll - sum(RTMB::dnorm(
          U[-1L, , drop = FALSE],
          phi * U[-n, , drop = FALSE],
          sd[1L] * innovation_sd,
          log = TRUE
        ))
      }
    }
  } else if (density_structure == "ou") {
    ## Match the state-space OU likelihood used by glmmTMB.cpp:603-608.
    ## This avoids dense dmvnorm() evaluation and preserves the Markov
    ## structure implied by adjacent time differences.
    nll <- 0
    time_diff <- diff(times)
    rho <- exp(-decay * time_diff)
    innovation_sd <- sd[1L] * sqrt(1 - rho * rho)

    if (simulation) {
      for (k in seq_len(reps)) {
        U_sim <- numeric(n)
        U_sim[1L] <- stats::rnorm(1L, mean = 0, sd = sd[1L])
        if (n > 1L) {
          for (j in 2:n) {
            U_sim[j] <- stats::rnorm(
              1L,
              mean = rho[j - 1L] * U_sim[j - 1L],
              sd = innovation_sd[j - 1L]
            )
          }
        }
        U_column <- U[, k]
        U_column[] <- U_sim
      }
    } else {
      nll <- -sum(RTMB::dnorm(U[1L, ], 0, sd[1L], log = TRUE))
      if (n > 1L) {
        nll <- nll - sum(RTMB::dnorm(
          U[-1L, , drop = FALSE],
          rho * U[-n, , drop = FALSE],
          innovation_sd,
          log = TRUE
        ))
      }
    }
  } else {
    ## Keep scale dimensions identical to t(U). A bare vector is ambiguous to
    ## RTMB::dmvnorm when there is exactly one block repetition.
    scale_matrix <- rep(sd, reps)
    dim(scale_matrix) <- c(n, reps)
    scale_matrix <- t(scale_matrix)
    nll <- -sum(RTMB::dmvnorm(t(U), Sigma = C, log = TRUE, scale = scale_matrix)
    )
  }

  ## Match C++ full-correlation reporting; equalto always reports its matrix.
  report_corr <- C
  if (name %in% c("ar1", "hetar1")) {
    report_corr <- matrix(phi, 1L, 1L)
  }
  if (name == "ou" && term$fullCor == 0) {
    report_corr <- matrix(decay, 1L, 1L)
  }
  if (name %in% c("exp", "gau", "mat") && term$fullCor == 0) {
    report_corr <- matrix(numeric(0), 0, 0)
  }

  conditional_full_cor <- c(
    "us", "cs", "homcs", "toep", "homtoep", "propto"
  )
  if (name %in% conditional_full_cor && term$fullCor == 0) {
    report_corr <- matrix(NaN, 1L, 1L)
  }

  report_sd <- if (name == "ar1" || (name == "ou" && term$fullCor == 0)) {
    sd[1L]
  } else {
    sd
  }
  list(
    nll = nll,
    corr = report_corr,
    sd = report_sd,
    fact_load = matrix(numeric(0), 0, 0),
    u = if (inherits(U, "simref")) NULL else as.vector(U)
  )
}
