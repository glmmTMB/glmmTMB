## Tests for RTMB spatial covariance structures

context("RTMB spatial covariance structures")

skip_if_not_installed("RTMB")

tol_logLik <- 1e-5
tol_fixef <- 1e-5
tol_varcorr <- 1e-4

set.seed(4001)
spatial_data <- expand.grid(
  x = c(0, 0.5, 2),
  y = c(0, 1),
  group = factor(seq_len(12))
)
spatial_data$pos <- glmmTMB::numFactor(spatial_data$x, spatial_data$y)
group_effect <- rnorm(nlevels(spatial_data$group), sd = 0.3)
linear_predictor <- 0.5 + 0.15 * spatial_data$x +
  group_effect[spatial_data$group]
spatial_data$gaussian_response <- linear_predictor +
  rnorm(nrow(spatial_data), sd = 0.5)
spatial_data$count <- rpois(nrow(spatial_data), exp(linear_predictor))

test_that("exponential spatial covariance matches TMB and distance formula", {
  spatial_sd <- 0.6
  spatial_range <- 1.4
  theta <- c(log(spatial_sd), log(spatial_range))
  theta_map <- factor(c(NA, NA))

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    gaussian_response ~ x + exp(0 + pos | group),
    family = gaussian,
    data = spatial_data,
    start = list(theta = theta),
    map = list(theta = theta_map),
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    gaussian_response ~ x + exp(0 + pos | group),
    family = gaussian,
    data = spatial_data,
    start = list(theta = theta),
    map = list(theta = theta_map),
    se = FALSE
  )

  correlation <- attr(VarCorr(m_rtmb)$cond$group, "correlation")
  coordinates <- glmmTMB::parseNumLevels(rownames(correlation))
  expected <- exp(-as.matrix(dist(coordinates)) / spatial_range)
  dimnames(expected) <- dimnames(correlation)

  expect_equal(
    as.numeric(logLik(m_rtmb)),
    as.numeric(logLik(m_tmb)),
    tolerance = tol_logLik
  )
  expect_equal(
    fixef(m_rtmb)$cond,
    fixef(m_tmb)$cond,
    tolerance = tol_fixef
  )
  expect_equal(
    VarCorr(m_rtmb),
    VarCorr(m_tmb),
    tolerance = tol_varcorr
  )
  expect_equal(correlation, expected, tolerance = tol_varcorr)
})

test_that("Gaussian spatial covariance matches TMB and distance formula", {
  spatial_sd <- 0.6
  spatial_range <- 1.4
  theta <- c(log(spatial_sd), log(spatial_range))
  theta_map <- factor(c(NA, NA))

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    gaussian_response ~ x + gau(0 + pos | group),
    family = gaussian,
    data = spatial_data,
    start = list(theta = theta),
    map = list(theta = theta_map),
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    gaussian_response ~ x + gau(0 + pos | group),
    family = gaussian,
    data = spatial_data,
    start = list(theta = theta),
    map = list(theta = theta_map),
    se = FALSE
  )

  correlation <- attr(VarCorr(m_rtmb)$cond$group, "correlation")
  coordinates <- glmmTMB::parseNumLevels(rownames(correlation))
  expected <- exp(
    -(as.matrix(dist(coordinates)) / spatial_range)^2
  )
  dimnames(expected) <- dimnames(correlation)

  expect_equal(
    as.numeric(logLik(m_rtmb)),
    as.numeric(logLik(m_tmb)),
    tolerance = tol_logLik
  )
  expect_equal(
    fixef(m_rtmb)$cond,
    fixef(m_tmb)$cond,
    tolerance = tol_fixef
  )
  expect_equal(
    VarCorr(m_rtmb),
    VarCorr(m_tmb),
    tolerance = tol_varcorr
  )
  expect_equal(correlation, expected, tolerance = tol_varcorr)
})

test_that("Matern spatial covariance matches TMB and distance formula", {
  spatial_sd <- 0.6
  spatial_range <- 1.4
  smoothness <- 1.2
  theta <- c(
    log(spatial_sd),
    log(spatial_range),
    log(smoothness)
  )
  theta_map <- factor(rep(NA, length(theta)))

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    gaussian_response ~ x + mat(0 + pos | group),
    family = gaussian,
    data = spatial_data,
    start = list(theta = theta),
    map = list(theta = theta_map),
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    gaussian_response ~ x + mat(0 + pos | group),
    family = gaussian,
    data = spatial_data,
    start = list(theta = theta),
    map = list(theta = theta_map),
    se = FALSE
  )

  correlation <- attr(VarCorr(m_rtmb)$cond$group, "correlation")
  coordinates <- glmmTMB::parseNumLevels(rownames(correlation))
  scaled_distance <- as.matrix(dist(coordinates)) / spatial_range
  expected <- matrix(1, nrow(scaled_distance), ncol(scaled_distance))
  off_diagonal <- row(expected) != col(expected)
  expected[off_diagonal] <-
    scaled_distance[off_diagonal]^smoothness *
      besselK(scaled_distance[off_diagonal], smoothness) /
      (exp(lgamma(smoothness)) * 2^(smoothness - 1))
  dimnames(expected) <- dimnames(correlation)

  expect_equal(
    as.numeric(logLik(m_rtmb)),
    as.numeric(logLik(m_tmb)),
    tolerance = tol_logLik
  )
  expect_equal(
    fixef(m_rtmb)$cond,
    fixef(m_tmb)$cond,
    tolerance = tol_fixef
  )
  expect_equal(
    VarCorr(m_rtmb),
    VarCorr(m_tmb),
    tolerance = tol_varcorr
  )
  expect_equal(correlation, expected, tolerance = tol_varcorr)

  set.seed(4004)
  simulation <- m_rtmb$obj$simulate(complete = TRUE)$yobs
  expect_length(simulation, nrow(spatial_data))
  expect_true(all(is.finite(simulation)))
})

test_that("exponential spatial covariance works with Poisson responses", {
  theta <- c(log(0.5), log(1.2))
  theta_map <- factor(c(NA, NA))

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x + exp(0 + pos | group),
    family = poisson,
    data = spatial_data,
    start = list(theta = theta),
    map = list(theta = theta_map),
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x + exp(0 + pos | group),
    family = poisson,
    data = spatial_data,
    start = list(theta = theta),
    map = list(theta = theta_map),
    se = FALSE
  )

  expect_equal(
    as.numeric(logLik(m_rtmb)),
    as.numeric(logLik(m_tmb)),
    tolerance = tol_logLik
  )
  expect_equal(
    fixef(m_rtmb)$cond,
    fixef(m_tmb)$cond,
    tolerance = tol_fixef
  )
  expect_equal(
    VarCorr(m_rtmb),
    VarCorr(m_tmb),
    tolerance = tol_varcorr
  )
})

test_that("Gaussian spatial covariance works in zero-inflation model", {
  zi_data <- spatial_data
  zi_data$count[seq(1, nrow(zi_data), by = 5)] <- 0
  thetazi <- c(log(0.5), log(1.2))
  theta_map <- factor(c(NA, NA))

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x,
    ziformula = ~ 1 + gau(0 + pos | group),
    family = poisson,
    data = zi_data,
    start = list(thetazi = thetazi),
    map = list(thetazi = theta_map),
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x,
    ziformula = ~ 1 + gau(0 + pos | group),
    family = poisson,
    data = zi_data,
    start = list(thetazi = thetazi),
    map = list(thetazi = theta_map),
    se = FALSE
  )

  expect_equal(
    as.numeric(logLik(m_rtmb)),
    as.numeric(logLik(m_tmb)),
    tolerance = tol_logLik
  )
  expect_equal(
    fixef(m_rtmb)$zi,
    fixef(m_tmb)$zi,
    tolerance = tol_fixef
  )
  expect_equal(
    VarCorr(m_rtmb),
    VarCorr(m_tmb),
    tolerance = tol_varcorr
  )
})

test_that("spatial covariance simulation works under RTMB", {
  theta <- c(log(0.5), log(1.2))
  theta_map <- factor(c(NA, NA))

  local_useRTMB(TRUE)
  exponential_model <- glmmTMB(
    gaussian_response ~ x + exp(0 + pos | group),
    family = gaussian,
    data = spatial_data,
    start = list(theta = theta),
    map = list(theta = theta_map),
    se = FALSE
  )
  gaussian_model <- glmmTMB(
    gaussian_response ~ x + gau(0 + pos | group),
    family = gaussian,
    data = spatial_data,
    start = list(theta = theta),
    map = list(theta = theta_map),
    se = FALSE
  )

  set.seed(4002)
  exp_simulation <- exponential_model$obj$simulate(complete = TRUE)$yobs
  set.seed(4003)
  gau_simulation <- gaussian_model$obj$simulate(complete = TRUE)$yobs

  expect_length(exp_simulation, nrow(spatial_data))
  expect_length(gau_simulation, nrow(spatial_data))
  expect_true(all(is.finite(exp_simulation)))
  expect_true(all(is.finite(gau_simulation)))
})

test_that("RTMB covariance structures validate parameter counts", {
  n <- 4L
  expected <- c(
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
    rr = 2L * n - 1L,
    propto = n * (n + 1L) / 2L + 1L,
    equalto = n * (n + 1L) / 2L
  )

  for (name in names(expected)) {
    expected_message <- if (name == "rr") {
      "Invalid covariance parameter count for 'rr'"
    } else {
      paste0(
        "Expected ", expected[[name]],
        " covariance parameters for '", name, "'"
      )
    }
    term <- list(
      blockCode = glmmTMB:::.valid_covstruct[name],
      blockSize = n,
      blockReps = 1L
    )

    expect_error(
      glmmTMB:::termwise_nll(
        numeric(n),
        numeric(expected[[name]] - 1L),
        term
      ),
      expected_message,
      fixed = TRUE
    )
  }
})
