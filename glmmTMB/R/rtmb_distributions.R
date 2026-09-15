## RTMB density and simulation helpers, including zero inflation and truncation.

vec_rtmb <- function(z, x) {
  rep(as.vector(z), length.out = length(x))
}

log_var_minus_mu_rtmb <- function(family_name, log_mu, etadisp, psi) {
  switch(
    family_name,
    nbinom1 = log_mu + etadisp,
    truncated_nbinom1 = log_mu + etadisp,
    nbinom2 = 2 * log_mu - etadisp,
    truncated_nbinom2 = 2 * log_mu - etadisp,
    nbinom12 = {
      log_mu_vec <- log_mu[seq_along(log_mu)]
      etadisp_vec <- etadisp[seq_along(etadisp)]
      log_mu_vec + RTMB::logspace_add(etadisp_vec, log_mu_vec - psi[1L])
    },
    stop("log(var - mu) not defined for distribution: ", family_name)
  )
}

#' Simulate from a zero-inflated density wrapper
#'
#' This helper implements the simulation branch used by [dZI()].  The
#' zero-inflation predictor `eta_zi` is on the logit scale, so the structural
#' zero probability is `p_zi = 1 / (1 + exp(-eta_zi))` and the complementary
#' conditional probability is `1 - p_zi = 1 / (1 + exp(eta_zi))`.  Simulated
#' values therefore come from the same mixture represented in the likelihood:
#' either the structural-zero component returns zero, or the wrapped
#' conditional density simulates the response.  Keeping this logic outside
#' [dZI()] makes the likelihood code easier to read and preserves RTMB's
#' convention that simulation is triggered by evaluating log-density functions
#' on `simref` objects.
#'
#' @param density A log-density function such as `RTMB::dpois()` or
#'   `RTMB::dnorm()` that also supports RTMB simulation through `simref`
#'   objects.
#' @param x The response vector, normally an RTMB `simref` object during
#'   simulation.
#' @param ... Distribution-specific arguments passed to `density`.
#' @param eta_zi Zero-inflation linear predictor on the logit scale.
#'
#' @return A zero vector on the log-density scale, after mutating the `simref`
#'   response object with simulated values.
#'
#' @noRd
simZI <- function(density, x, ..., eta_zi) {
  prob_nonzero <- 1 / (1 + exp(eta_zi))

  nonzero <- as.logical(stats::rbinom(length(x), 1, prob_nonzero))
  density_args <- list(...)
  density_args <- lapply(density_args, function(arg) {
    if (length(arg) == length(x)) arg[nonzero] else arg
  })

  if (any(nonzero)) {
    do.call(
      density,
      c(list(x = x[nonzero]), density_args, list(log = TRUE))
    )
  }
  if (any(!nonzero)) {
    structural_zero <- x[!nonzero]
    structural_zero[] <- 0
  }
  rep(0, length(x))
}

#' Add zero inflation to a density function
#'
#' `dZI()` takes an ordinary density function and
#' returns a new density function with glmmTMB-style zero-inflation behavior.
#' When `eta_zi` is `NULL`, the wrapper deliberately reduces to the original
#' density, so non-zero-inflated models use the same likelihood path.  When
#' `eta_zi` is supplied, the likelihood is the standard zero-inflated mixture.
#' A nonzero observation can only come from the conditional distribution, so it
#' contributes `log(1 - p_zi) + log f(y)`.  An observed zero can come either
#' from the structural-zero component or from the conditional distribution, so
#' it contributes `log(p_zi + (1 - p_zi) * f(0))`.  The wrapper evaluates these
#' terms on the log scale, using `RTMB::logspace_add()` for the zero case so
#' the two possible zero sources are combined stably.  The separate `is_zero`
#' argument tells the wrapper which observations should use the zero-mixture
#' formula; this is especially important for truncated or hurdle-like families
#' where the conditional density should not be evaluated at zero.
#'
#' @param density A log-density function to wrap.
#'
#' @return A function with the same distribution-specific arguments as
#'   `density`, plus `eta_zi`, `log`, and `is_zero`.
#'
#' @noRd
dZI <- function(density) {
  force(density)

  function(x, ..., eta_zi = NULL, log = FALSE, is_zero = NULL) {
    if (inherits(x, "simref") && !is.null(eta_zi)) {
      return(simZI(density, x, ..., eta_zi = eta_zi))
    }

    x <- osa_value(x)
    if (is.null(eta_zi)) {
      loglik <- density(x, ..., log = TRUE)
      if (log) {
        return(loglik)
      }
      return(exp(loglik))
    }

    if (is.null(is_zero)) {
      is_zero <- x == 0
    }
    has_zero <- any(is_zero)

    loglik <- density(x, ..., log = TRUE)
    log_1mpz <- -RTMB::logspace_add(0, eta_zi)
    ans <- log_1mpz + loglik

    if (has_zero) {
      log_pz <- -RTMB::logspace_add(0, -eta_zi)
      ans[is_zero] <- RTMB::logspace_add(
        log_pz[is_zero],
        ans[is_zero]
      )
    }
    if (log) ans else exp(ans)
  }
}

## Fitting uses RTMB::dbinom_robust() to match glmmTMB.cpp:981.
## Simulation follows glmmTMB.cpp:982 but delegates to RTMB::dbinom(),
## because RTMB::dbinom_robust() does not provide rbinom_robust().
dbinom_robust_rtmb <- function(x, size, logit_p, log = FALSE) {
  if (inherits(x, "simref")) {
    prob <- 1 / (1 + exp(-logit_p))
    return(RTMB::dbinom(x, size = size, prob = prob, log = log))
  }
  RTMB::dbinom_robust(x, size = size, logit_p = logit_p, log = log)
}

## Translated from glmmtmb::logspace_gamma(), distrib.h:27-47.
## This avoids lgamma(exp(x)) underflow for very small x, matching the C++
## stable beta-binomial density helper.
logspace_gamma_rtmb <- RTMB::Vectorize(
  function(x) {
    `if` <- RTMB::ADoverload("if")
    if (x < -150) -x else lgamma(exp(x))
  },
  vectorize.args = "x"
)

## Translated from glmmtmb::dbetabinom_robust(), distrib.h:49-62, and
## the betabinomial_family case in glmmTMB.cpp:1037-1047.
dbetabinom_robust_rtmb <- function(x, log_shape1, log_shape2, size,
                                   log = FALSE) {
  if (inherits(x, "simref")) {
    log_shape1 <- vec_rtmb(log_shape1, x)
    log_shape2 <- vec_rtmb(log_shape2, x)
    size <- vec_rtmb(size, x)
    ans <- numeric(length(x))
    for (i in seq_along(ans)) {
      prob <- stats::rbeta(1L, exp(log_shape1[i]), exp(log_shape2[i]))
      ans[i] <- stats::rbinom(1L, size = size[i], prob = prob)
    }
    x[] <- ans
    return(rep(0, length(x)))
  }

  log_x <- log(x)
  log_size_minus_x <- log(size - x)
  ans <-
    lgamma(size + 1) - lgamma(x + 1) - lgamma(size - x + 1) +
    logspace_gamma_rtmb(RTMB::logspace_add(log_x, log_shape1)) +
    logspace_gamma_rtmb(RTMB::logspace_add(log_size_minus_x, log_shape2)) -
    lgamma(size + exp(log_shape1) + exp(log_shape2)) +
    lgamma(exp(log_shape1) + exp(log_shape2)) -
    logspace_gamma_rtmb(log_shape1) - logspace_gamma_rtmb(log_shape2)
  if (log) ans else exp(ans)
}

## Mean-parameterized Conway-Maxwell-Binomial density; translated from
## dcombinom2() and combinom_utils in TMB's distributions_R.hpp:660-695 and
## tiny_ad/compois/combinom.hpp:1-96.
combinom_logZ_rtmb <- function(logit_p, nu, size) {
  ans <- -Inf
  for (k in 0:size) {
    ans <- RTMB::logspace_add(
      ans,
      nu * lchoose(size, k) + k * logit_p
    )
  }
  ans
}

combinom_moments_rtmb <- function(logit_p, nu, size) {
  log_z <- combinom_logZ_rtmb(logit_p, nu, size)
  mean <- second_moment <- 0
  for (k in 0:size) {
    probability <- exp(
      nu * lchoose(size, k) + k * logit_p - log_z
    )
    mean <- mean + k * probability
    second_moment <- second_moment + k * k * probability
  }
  list(mean = mean, variance = second_moment - mean * mean)
}

combinom_logitp_rtmb <- function(mean, nu, size) {
  "if" <- RTMB::ADoverload("if")
  lower <- -30
  upper <- 30
  for (i in seq_len(10L)) {
    midpoint <- (lower + upper) / 2
    midpoint_mean <- combinom_moments_rtmb(midpoint, nu, size)$mean
    new_lower <- if (midpoint_mean > mean) lower else midpoint
    new_upper <- if (midpoint_mean > mean) midpoint else upper
    lower <- new_lower
    upper <- new_upper
  }

  logit_p <- (lower + upper) / 2
  for (i in seq_len(5L)) {
    moments <- combinom_moments_rtmb(logit_p, nu, size)
    logit_p <- logit_p - (moments$mean - mean) / moments$variance
  }
  logit_p
}

rcombinom2_rtmb <- function(mean, nu, size) {
  logit_p <- combinom_logitp_rtmb(mean, nu, size)
  log_weights <- vapply(
    0:size,
    function(k) nu * lchoose(size, k) + k * logit_p,
    numeric(1)
  )
  weights <- exp(log_weights - max(log_weights))
  random_number <- stats::runif(1L) * sum(weights)
  which(cumsum(weights) >= random_number)[1L] - 1L
}

dcombinom2_rtmb <- function(x, size, mean, nu, log = FALSE) {
  if (inherits(x, "simref")) {
    mean <- as.vector(mean)
    nu <- as.vector(nu)
    ans <- numeric(length(x))
    for (i in seq_along(ans)) {
      ans[i] <- rcombinom2_rtmb(mean[i], nu[i], size[i])
    }
    x[] <- ans
    return(rep(0, length(x)))
  }

  "[<-" <- RTMB::ADoverload("[<-")
  ans <- mean * 0
  for (i in seq_along(x)) {
    logit_p <- combinom_logitp_rtmb(mean[i], nu[i], size[i])
    ans[i] <- nu[i] * lchoose(size[i], x[i]) + x[i] * logit_p -
      combinom_logZ_rtmb(logit_p, nu[i], size[i])
  }
  if (log) ans else exp(ans)
}

## Translated from the Gamma_family case in glmmTMB.cpp:991-996 and
## zt_lik_zero(), glmmTMB.cpp:959. Gamma is mean/shape parameterized here:
## shape = phi and scale = mu / phi. Exact zeros are excluded from the
## conditional density so zero-inflated Gamma behaves as a hurdle model.
dgamma_rtmb <- function(x, mean, shape, log = FALSE) {
  if (inherits(x, "simref")) {
    x[] <- stats::rgamma(
      length(x),
      shape = as.vector(shape),
      scale = as.vector(mean / shape)
    )
    return(rep(0, length(x)))
  }

  is_zero <- x == 0
  if (!any(is_zero)) {
    ans <- RTMB::dgamma(x, shape = shape, scale = mean / shape, log = TRUE)
  } else {
    not_zero <- !is_zero
    if (!any(not_zero)) {
      return(rep(if (log) -Inf else 0, length(x)))
    }

    "[<-" <- RTMB::ADoverload("[<-")
    subset_arg <- function(arg) {
      if (length(arg) == length(x)) arg[not_zero] else arg
    }
    mean_nz <- subset_arg(mean)
    shape_nz <- subset_arg(shape)
    loglik_nz <- RTMB::dgamma(
      x[not_zero],
      shape = shape_nz,
      scale = mean_nz / shape_nz,
      log = TRUE
    )
    ans <- rep(loglik_nz[1L] * 0 - Inf, length(x))
    ans[not_zero] <- loglik_nz
  }
  if (log) ans else exp(ans)
}

## Translated from the beta_family case in glmmTMB.cpp:997-1002 and
## zt_lik_zero(), glmmTMB.cpp:959. The beta family uses Ferrari-Cribari-Neto
## parameterization: shape1 = mu * phi and shape2 = (1 - mu) * phi. Exact zeros
## are excluded from the conditional density so zero-inflated beta behaves as a
## hurdle model.
dbeta_rtmb <- function(x, mean, phi, log = FALSE) {
  shape1 <- mean * phi
  shape2 <- (1 - mean) * phi

  if (inherits(x, "simref")) {
    x[] <- stats::rbeta(
      length(x),
      shape1 = as.vector(shape1),
      shape2 = as.vector(shape2)
    )
    return(rep(0, length(x)))
  }

  is_zero <- x == 0
  if (!any(is_zero)) {
    ans <- RTMB::dbeta(x, shape1 = shape1, shape2 = shape2, log = TRUE)
  } else {
    not_zero <- !is_zero
    if (!any(not_zero)) {
      return(rep(if (log) -Inf else 0, length(x)))
    }

    "[<-" <- RTMB::ADoverload("[<-")
    subset_arg <- function(arg) {
      if (length(arg) == length(x)) arg[not_zero] else arg
    }
    shape1_nz <- subset_arg(shape1)
    shape2_nz <- subset_arg(shape2)
    loglik_nz <- RTMB::dbeta(
      x[not_zero],
      shape1 = shape1_nz,
      shape2 = shape2_nz,
      log = TRUE
    )
    ans <- rep(loglik_nz[1L] * 0 - Inf, length(x))
    ans[not_zero] <- loglik_nz
  }
  if (log) ans else exp(ans)
}

## Translated from the ordbeta_family case in glmmTMB.cpp:1004-1031.
dordbeta_rtmb <- function(x, eta, mean, phi, cutpoints, log = FALSE) {
  if (inherits(x, "simref")) {
    eta <- as.vector(eta)
    mean <- as.vector(mean)
    phi <- as.vector(phi)
    ans <- numeric(length(x))
    for (i in seq_along(ans)) {
      if (stats::runif(1L) < stats::plogis(cutpoints[1L] - eta[i])) {
        ans[i] <- 0
      } else if (stats::runif(1L) < stats::plogis(eta[i] - cutpoints[2L])) {
        ans[i] <- 1
      } else {
        ans[i] <- stats::rbeta(
          1L,
          shape1 = mean[i] * phi[i],
          shape2 = (1 - mean[i]) * phi[i]
        )
      }
    }
    x[] <- ans
    return(rep(0, length(x)))
  }

  "[<-" <- RTMB::ADoverload("[<-")
  is_zero <- x == 0
  is_one <- x == 1
  is_interior <- !is_zero & !is_one
  ans <- eta * 0

  if (any(is_zero)) {
    ans[is_zero] <- -RTMB::logspace_add(
      0,
      eta[is_zero] - cutpoints[1L]
    )
  }
  if (any(is_one)) {
    ans[is_one] <- -RTMB::logspace_add(
      0,
      cutpoints[2L] - eta[is_one]
    )
  }
  if (any(is_interior)) {
    log_lower <- -RTMB::logspace_add(
      0,
      cutpoints[1L] - eta[is_interior]
    )
    log_upper <- -RTMB::logspace_add(
      0,
      cutpoints[2L] - eta[is_interior]
    )
    log_middle <- RTMB::logspace_sub(log_lower, log_upper)
    ans[is_interior] <- log_middle + RTMB::dbeta(
      x[is_interior],
      shape1 = mean[is_interior] * phi[is_interior],
      shape2 = (1 - mean[is_interior]) * phi[is_interior],
      log = TRUE
    )
  }
  if (log) ans else exp(ans)
}

## Translated from the lognormal_family case in glmmTMB.cpp:1164-1179.
## The lognormal family is parameterized by mean and SD on the data scale.
dlognormal_rtmb <- function(x, mean, sd, log = FALSE) {
  log_var <- RTMB::logspace_add(2 * (log(sd) - log(mean)), 0)
  meanlog <- log(mean) - log_var / 2
  sdlog <- sqrt(log_var)

  if (inherits(x, "simref")) {
    x[] <- stats::rlnorm(length(x), meanlog = as.vector(meanlog),
                         sdlog = as.vector(sdlog))
    return(rep(0, length(x)))
  }

  is_zero <- x == 0
  if (!any(is_zero)) {
    ans <- RTMB::dnorm(log(x), meanlog, sdlog, log = TRUE) - log(x)
  } else {
    not_zero <- !is_zero
    if (!any(not_zero)) {
      return(rep(if (log) -Inf else 0, length(x)))
    }

    "[<-" <- RTMB::ADoverload("[<-")
    subset_arg <- function(arg) {
      if (length(arg) == length(x)) arg[not_zero] else arg
    }
    meanlog_nz <- subset_arg(meanlog)
    sdlog_nz <- subset_arg(sdlog)
    loglik_nz <- RTMB::dnorm(
      log(x[not_zero]),
      meanlog_nz,
      sdlog_nz,
      log = TRUE
    ) - log(x[not_zero])
    ans <- rep(loglik_nz[1L] * 0 - Inf, length(x))
    ans[not_zero] <- loglik_nz
  }
  if (log) ans else exp(ans)
}

## Translated from the t_family case in glmmTMB.cpp:1182-1190.
## The response is standardized by the fitted scale phi, so the log-density
## subtracts log(phi), represented by etadisp in the C++ code.
dt_rtmb <- function(x, mean, scale, df, log = FALSE) {
  if (inherits(x, "simref")) {
    x[] <- as.vector(mean) + as.vector(scale) * stats::rt(length(x), df)
    return(rep(0, length(x)))
  }

  ans <- RTMB::dt((x - mean) / scale, df = df, log = TRUE) - log(scale)
  if (log) ans else exp(ans)
}

## Translated from glmmtmb::dskewnorm(), distrib.h:78-88, and the
## skewnormal_family case in glmmTMB.cpp:975-980.
dskewnormal_rtmb <- function(x, mean, sd, alpha, log = FALSE) {
  delta <- alpha / sqrt(1 + alpha^2)
  omega <- sd / sqrt(1 - 2 / pi * delta^2)
  xi <- mean - omega * delta * sqrt(2 / pi)

  if (inherits(x, "simref")) {
    n <- length(x)
    xi <- vec_rtmb(xi, x)
    omega <- vec_rtmb(omega, x)
    delta <- vec_rtmb(delta, x)
    ans <- numeric(n)
    for (i in seq_len(n)) {
      chi <- abs(stats::rnorm(1L))
      nrv <- stats::rnorm(1L)
      z <- delta[i] * chi + sqrt(1 - delta[i]^2) * nrv
      ans[i] <- xi[i] + omega[i] * z
    }
    x[] <- ans
    return(rep(0, length(x)))
  }

  z <- (x - xi) / omega
  ans <- log(2) - log(omega) +
    RTMB::dnorm(z, 0, 1, log = TRUE) +
    log(RTMB::pnorm(alpha * z))
  if (log) ans else exp(ans)
}

## Fitting translates glmmTMB.cpp:1042-1075 for nbinom1/nbinom2:
## both families use dnbinom_robust(log_mu, log_var_minus_mu). Simulation
## follows the same mean/variance by converting back to size/mu.
dnbinom_robust_rtmb <- function(x, log_mu, log_var_minus_mu, log = FALSE) {
  if (inherits(x, "simref")) {
    mu <- exp(as.vector(log_mu))
    size <- exp(as.vector(2 * log_mu - log_var_minus_mu))
    x[] <- stats::rnbinom(length(x), size = size, mu = mu)
    return(rep(0, length(x)))
  }
  RTMB::dnbinom_robust(x, log_mu = log_mu, log_var_minus_mu = log_var_minus_mu,
                       log = log)
}

## Translated from glmmtmb::dgenpois(); distrib.h:64-78.
dgenpois_log_rtmb <- RTMB::Vectorize(
  function(x, theta, lambda) {
    term <- theta + lambda * x
    log(theta) + (x - 1) * log(term) - term - lgamma(x + 1)
  },
  vectorize.args = c("x", "theta", "lambda")
)

## Translated from glmmtmb::rgenpois(); distrib.h:172-183.
rgenpois_rtmb <- function(theta, lambda) {
  ans <- 0
  random_number <- stats::runif(1L)
  kum <- exp(dgenpois_log_rtmb(0, theta, lambda))
  while (random_number > kum) {
    ans <- ans + 1
    kum <- kum + exp(dgenpois_log_rtmb(ans, theta, lambda))
  }
  ans
}

## Translated from glmmtmb::rtruncated_genpois(); distrib.h:186-197.
rtruncated_genpois_rtmb <- function(theta, lambda) {
  nloop <- 10000L
  counter <- 0L
  ans <- rgenpois_rtmb(theta, lambda)
  while (ans < 1 && counter < nloop) {
    ans <- rgenpois_rtmb(theta, lambda)
    counter <- counter + 1L
  }
  if (counter == nloop && ans < 1) {
    warning(
      "Simulation of zero-truncated data reached the maximum number of retries ",
      "and still returned zeros. ",
      "Possibly due to low estimated mean.",
      call. = FALSE
    )
  }
  ans
}

## Translated from the genpois_family case in glmmTMB.cpp:1128-1133.
dgenpois_rtmb <- function(x, theta, lambda, log = FALSE) {
  if (inherits(x, "simref")) {
    theta <- vec_rtmb(theta, x)
    lambda <- vec_rtmb(lambda, x)
    ans <- numeric(length(x))
    for (i in seq_along(ans)) {
      ans[i] <- rgenpois_rtmb(theta[i], lambda[i])
    }
    x[] <- ans
    return(rep(0, length(x)))
  }

  ans <- dgenpois_log_rtmb(x, theta, lambda)
  if (log) ans else exp(ans)
}

## AD-compatible Lambert W transformation used by the Bell family; translated
## from glmmtmb::LambertW(), distrib.h:486-521.
lambertW_rtmb <- RTMB::Vectorize(
  function(x) {
    "if" <- RTMB::ADoverload("if")
    logx <- log(x)
    y <- if (logx > 0) logx else x * 0
    for (i in seq_len(100L)) {
      y <- y - (y - exp(logx - y)) / (1 + y)
    }
    y
  },
  vectorize.args = "x"
)

## Translation of glmmtmb::Bell(), distrib.h:442-464.
bell_number_rtmb <- function(n) {
  if (n < 2L) {
    return(1)
  }

  bell <- bell_new <- numeric(n)
  bell[1L] <- 1
  for (i in seq_len(n - 1L)) {
    bell_new[1L] <- bell[i]
    for (j in seq_len(i)) {
      bell_new[j + 1L] <- bell[j] + bell_new[j]
    }
    bell <- bell_new
  }
  bell_new[n]
}

## Translation of glmmtmb::rbell() and glmmtmb::dbell(),
## distrib.h:417-432 and 466-478.
dbell_rtmb <- function(x, mean, log = FALSE) {
  theta <- lambertW_rtmb(mean)
  if (inherits(x, "simref")) {
    theta <- as.vector(theta)
    ans <- numeric(length(x))
    for (i in seq_along(ans)) {
      n_compound <- stats::rpois(1L, expm1(theta[i]))
      if (n_compound > 0L) {
        ans[i] <- sum(vapply(
          seq_len(n_compound),
          function(j) rtruncated_poisson_rtmb(theta[i]),
          numeric(1)
        ))
      }
    }
    x[] <- ans
    return(rep(0, length(x)))
  }

  bell_number <- vapply(as.integer(x), bell_number_rtmb, numeric(1))
  ans <- x * log(theta) - exp(theta) + 1 +
    log(bell_number) - lgamma(x + 1)
  if (log) ans else exp(ans)
}

## Translated from calc_log_nzprob(), glmmTMB.cpp:286-291.
log_nzprob_truncated_genpois_rtmb <- RTMB::Vectorize(
  function(theta) RTMB::logspace_sub(0, -theta),
  vectorize.args = "theta"
)

## Translated from the truncated_genpois_family case,
## glmmTMB.cpp:1134-1139.
dtruncated_genpois_rtmb <- function(x, theta, lambda, log = FALSE) {
  if (inherits(x, "simref")) {
    theta <- vec_rtmb(theta, x)
    lambda <- vec_rtmb(lambda, x)
    ans <- numeric(length(x))
    for (i in seq_along(ans)) {
      ans[i] <- rtruncated_genpois_rtmb(theta[i], lambda[i])
    }
    x[] <- ans
    return(rep(0, length(x)))
  }

  log_nzprob <- log_nzprob_truncated_genpois_rtmb(theta)
  ans <- dgenpois_log_rtmb(x, theta, lambda) - log_nzprob

  is_zero <- x < 0.001
  if (any(is_zero)) {
    ans[is_zero] <- -Inf
  }
  if (log) ans else exp(ans)
}

## Translated from the compois_family case in glmmTMB.cpp:1115-1119.
## Fitting uses RTMB::dcompois2(mean, nu); simulation uses RTMB's internal
## rcompois2() because the exported density is the fitting interface.
dcompois2_rtmb <- function(x, mean, nu, log = FALSE) {
  if (inherits(x, "simref")) {
    x[] <- get("rcompois2", envir = asNamespace("RTMB"))(
      length(x),
      mean = as.vector(mean),
      nu = as.vector(nu)
    )
    return(rep(0, length(x)))
  }
  RTMB::dcompois2(x, mean = mean, nu = nu, log = log)
}

## Translation of glmmtmb::rtruncated_compois2(); distrib.h:280-289.
rtruncated_compois2_rtmb <- function(n, mean, nu) {
  rcompois2 <- get("rcompois2", envir = asNamespace("RTMB"))
  mean <- rep(mean, length.out = n)
  nu <- rep(nu, length.out = n)
  ans <- rcompois2(n, mean = mean, nu = nu)

  nloop <- 10000L
  counter <- 0L
  while (any(ans < 1) && counter < nloop) {
    zero <- ans < 1
    ans[zero] <- rcompois2(sum(zero), mean = mean[zero], nu = nu[zero])
    counter <- counter + 1L
  }
  if (counter == nloop && any(ans < 1)) {
    warning(
      "Simulation of zero-truncated data reached the maximum number of retries ",
      "and still returned zeros. ",
      "Possibly due to low estimated mean.",
      call. = FALSE
    )
  }
  ans
}

## translation of glmmtmb::rtruncated_nbinom(); distrib.h:130-168
rtruncated_nbinom_rtmb <- function(n, size, k = 0L, mu) {
  ans <- numeric(n)
  size <- rep(size, length.out = n)
  mu <- rep(mu, length.out = n)

  if (any(size <= 0)) {
    stop("non-positive size in k-truncated-neg-bin simulator")
  }
  if (any(mu <= 0)) {
    stop("non-positive mu in k-truncated-neg-bin simulator")
  }
  if (k < 0) {
    stop("negative k in k-truncated-neg-bin simulator")
  }

  for (i in seq_len(n)) {
    p <- size[i] / (mu[i] + size[i])
    q <- mu[i] / (mu[i] + size[i])
    m <- ceiling(max((k + 1) * p - size[i] * q, 0))

    repeat {
      x <- stats::rnbinom(1L, size = size[i] + m, prob = p) + m
      if (m > 0) {
        a <- 1
        u <- stats::runif(1L)
        for (j in seq_len(m)) {
          a <- a * (k + 2 - j) / (x - j + 1)
        }
        if (u < a && x > k) {
          break
        }
      } else if (x > k) {
        break
      }
    }
    ans[i] <- x
  }
  ans
}

## Scalar formulas vectorized to mirror the observation loop in
## calc_log_nzprob(), glmmTMB.cpp:270-290
log_nzprob_truncated_poisson_rtmb <- RTMB::Vectorize(
  function(mu) RTMB::logspace_sub(0, -mu),
  vectorize.args = "mu"
)

log_nzprob_truncated_nbinom1_rtmb <- RTMB::Vectorize(
  function(mu, log_phi) {
    log_phi_plus_one <- RTMB::logspace_add(0, log_phi)
    RTMB::logspace_sub(0, -mu / exp(log_phi) * log_phi_plus_one)
  },
  vectorize.args = c("mu", "log_phi")
)

log_nzprob_truncated_nbinom2_rtmb <- RTMB::Vectorize(
  function(log_mu, log_size) {
    log_ratio_plus_one <- RTMB::logspace_add(0, log_mu - log_size)
    RTMB::logspace_sub(0, -exp(log_size) * log_ratio_plus_one)
  },
  vectorize.args = c("log_mu", "log_size")
)

log_nzprob_truncated_compois_rtmb <- RTMB::Vectorize(
  function(mean, nu) {
    RTMB::logspace_sub(
      0,
      RTMB::dcompois2(0, mean = mean, nu = nu, log = TRUE)
    )
  },
  vectorize.args = c("mean", "nu")
)

## Translation of glmmtmb::rtruncated_poisson(); distrib.h:94-128.
rtruncated_poisson_rtmb <- function(mu, k = 0L) {
  if (mu <= 0) {
    stop("non-positive mu in k-truncated-poisson simulator")
  }
  if (k < 0) {
    stop("negative k in k-truncated-poisson simulator")
  }

  mdoub <- max(k + 1 - mu, 0)
  m <- ceiling(mdoub)

  repeat {
    x <- stats::rpois(1L, mu) + m
    if (m > 0) {
      a <- 1
      u <- stats::runif(1L)
      for (j in seq_len(m) - 1L) {
        a <- a * (k + 1 - j) / (x - j)
      }
      if (u < a && x > k) {
        return(x)
      }
    } else if (x > k) {
      return(x)
    }
  }
}

## Translated from calc_log_nzprob() and truncated_compois_family,
## glmmTMB.cpp:293-294 and 1121-1127.
dtruncated_compois2_rtmb <- function(x, mean, nu, log = FALSE) {
  if (inherits(x, "simref")) {
    x[] <- rtruncated_compois2_rtmb(
      length(x),
      mean = as.vector(mean),
      nu = as.vector(nu)
    )
    return(rep(0, length(x)))
  }

  log_nzprob <- log_nzprob_truncated_compois_rtmb(mean, nu)
  ans <- RTMB::dcompois2(x, mean = mean, nu = nu, log = TRUE) - log_nzprob

  is_zero <- x < 0.001
  if (any(is_zero)) {
    ans[is_zero] <- -Inf
  }
  if (log) ans else exp(ans)
}

## Translated from calc_log_nzprob() and the truncated_nbinom1_family
## likelihood case, glmmTMB.cpp:274-277 and 1042-1064
dtruncated_nbinom1_rtmb <- function(x, log_mu, log_var_minus_mu, log_phi,
                                    log = FALSE) {
  if (inherits(x, "simref")) {
    sim_mu <- exp(as.vector(log_mu))
    sim_phi <- exp(as.vector(log_phi))
    x[] <- rtruncated_nbinom_rtmb(
      length(x),
      size = sim_mu / sim_phi,
      k = 0L,
      mu = sim_mu
    )
    return(rep(0, length(x)))
  }

  mu <- exp(log_mu)
  log_nzprob <- log_nzprob_truncated_nbinom1_rtmb(mu, log_phi)
  ans <- RTMB::dnbinom_robust(x, log_mu = log_mu,
                              log_var_minus_mu = log_var_minus_mu,
                              log = TRUE) - log_nzprob

  is_zero <- x < 0.001
  if (any(is_zero)) {
    ans[is_zero] <- -Inf
  }
  if (log) ans else exp(ans)
}

## Translated from calc_log_nzprob() and the truncated_nbinom2_family
## likelihood case, glmmTMB.cpp:278-283 and 1066-1081
dtruncated_nbinom2_rtmb <- function(x, log_mu, log_var_minus_mu, log_size,
                                    log = FALSE) {
  if (inherits(x, "simref")) {
    x[] <- rtruncated_nbinom_rtmb(
      length(x),
      size = exp(as.vector(log_size)),
      k = 0L,
      mu = exp(as.vector(log_mu))
    )
    return(rep(0, length(x)))
  }

  log_nzprob <- log_nzprob_truncated_nbinom2_rtmb(log_mu, log_size)
  ans <- RTMB::dnbinom_robust(x, log_mu = log_mu,
                              log_var_minus_mu = log_var_minus_mu,
                              log = TRUE) - log_nzprob

  is_zero <- x < 0.001
  if (any(is_zero)) {
    ans[is_zero] <- -Inf
  }
  if (log) ans else exp(ans)
}

## zero-truncated poisson density
dtruncated_poisson_rtmb <- function(x, lambda, log = FALSE) {
  if (inherits(x, "simref")) {
    lambda <- vec_rtmb(lambda, x)
    ans <- numeric(length(x))
    for (i in seq_along(ans)) {
      ans[i] <- rtruncated_poisson_rtmb(lambda[i], k = 0L)
    }
    x[] <- ans
    return(rep(0, length(x)))
  }

  log_nzprob <- RTMB::logspace_sub(0, -lambda)
  ans <- RTMB::dpois(x, lambda = lambda, log = TRUE) - log_nzprob

  ## the conditional distribution has strictly positive support
  is_zero <- x < 0.001
  if (any(is_zero)) {
    ## return -Inf to let dZI() treat observed zeros as structural
    ans[is_zero] <- -Inf
  }
  if (log) ans else exp(ans)
}

dcauchy_rtmb <- function(x, location, scale, log = FALSE) {
  resid <- (x - location) / scale
  ans <- -log(pi) - log(scale) - log1p(resid * resid)
  if (log) ans else exp(ans)
}

dlkj_rtmb <- function(x, eta, log = FALSE) {
  "[<-" <- RTMB::ADoverload("[<-")

  len <- length(x)
  if (len == 0) {
    return(if (log) 0 else 1)
  }

  n <- (1 + sqrt(1 + 8 * len)) / 2
  if (abs(n - round(n)) > sqrt(.Machine$double.eps)) {
    stop("Invalid number of LKJ correlation parameters: ", len)
  }
  n <- as.integer(round(n))

  L <- diag(n)
  k <- 1L
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i > j) {
        L[i, j] <- x[k]
        k <- k + 1L
      }
    }
  }

  row_sums <- L * L
  log_det_x <- 0
  for (i in seq_len(n)) {
    log_det_x <- log_det_x - log(sum(row_sums[i, ]))
  }
  ans <- (eta - 1) * log_det_x
  if (log) ans else exp(ans)
}
