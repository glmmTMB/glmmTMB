## Test cases for the RTMB Tweedie family
## Fit each model with RTMB and legacy TMB, then compare likelihoods,
## parameters, predictions, and covariance estimates where applicable.

context("RTMB Tweedie backend")

skip_if_not_installed("RTMB")

tol_logLik <- 1e-5
tol_fixef <- 1e-5
tol_varcorr <- 1e-4

rtweedie_test <- function(n, mu, phi, p) {
  lambda <- mu^(2 - p) / (phi * (2 - p))
  alpha <- (2 - p) / (1 - p)
  scale <- phi * (p - 1) * mu^(p - 1)
  count <- rpois(n, lambda)
  rgamma(n, shape = -count * alpha, scale = scale)
}

set.seed(2101)
tweedie_dat <- expand.grid(
  obs = seq_len(5),
  g = factor(seq_len(30))
)
tweedie_dat$x <- rnorm(nrow(tweedie_dat))
tweedie_dat$z <- rnorm(nrow(tweedie_dat))
tweedie_dat$off <- runif(nrow(tweedie_dat), -0.15, 0.15)
tweedie_dat$w <- runif(nrow(tweedie_dat), 0.5, 1.5)
tweedie_effect <- rnorm(nlevels(tweedie_dat$g), sd = 0.2)
tweedie_mu <- exp(0.5 + 0.25 * tweedie_dat$x + tweedie_effect[tweedie_dat$g])
tweedie_phi <- 0.7
tweedie_power <- 1.5
tweedie_dat$y <- rtweedie_test(
  nrow(tweedie_dat),
  mu = tweedie_mu,
  phi = tweedie_phi,
  p = tweedie_power
)

tweedie_start <- list(psi = qlogis(tweedie_power - 1))
tweedie_map <- list(psi = factor(NA))

expect_tweedie_matches <- function(m_rtmb, m_tmb, check_varcorr = FALSE) {
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
  if (check_varcorr) {
    expect_equal(
      as.numeric(VarCorr(m_rtmb)$cond$g),
      as.numeric(VarCorr(m_tmb)$cond$g),
      tolerance = tol_varcorr
    )
  }
}

test_that("tweedie: fixed conditional effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    family = tweedie(),
    data = tweedie_dat,
    start = tweedie_start,
    map = tweedie_map,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    family = tweedie(),
    data = tweedie_dat,
    start = tweedie_start,
    map = tweedie_map,
    se = FALSE
  )

  expect_tweedie_matches(m_rtmb, m_tmb)
})

test_that("tweedie: offsets, weights, and dispersion effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x + offset(off),
    dispformula = ~ z,
    weights = w,
    family = tweedie(),
    data = tweedie_dat,
    start = tweedie_start,
    map = tweedie_map,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x + offset(off),
    dispformula = ~ z,
    weights = w,
    family = tweedie(),
    data = tweedie_dat,
    start = tweedie_start,
    map = tweedie_map,
    se = FALSE
  )

  expect_tweedie_matches(m_rtmb, m_tmb)
  expect_equal(fixef(m_rtmb)$disp, fixef(m_tmb)$disp, tolerance = tol_fixef)
})

test_that("tweedie: conditional random intercept", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x + (1 | g),
    family = tweedie(),
    data = tweedie_dat,
    start = tweedie_start,
    map = tweedie_map,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x + (1 | g),
    family = tweedie(),
    data = tweedie_dat,
    start = tweedie_start,
    map = tweedie_map,
    se = FALSE
  )

  expect_tweedie_matches(m_rtmb, m_tmb, check_varcorr = TRUE)
})

test_that("tweedie: estimated power parameter", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    family = tweedie(),
    data = tweedie_dat,
    start = tweedie_start,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    family = tweedie(),
    data = tweedie_dat,
    start = tweedie_start,
    se = FALSE
  )

  expect_tweedie_matches(m_rtmb, m_tmb)
  expect_equal(
    m_rtmb$obj$env$last.par.best["psi"],
    m_tmb$obj$env$last.par.best["psi"],
    tolerance = tol_fixef
  )
  expect_equal(
    family_params(m_rtmb),
    family_params(m_tmb),
    tolerance = tol_fixef
  )
})

test_that("tweedie: zero inflation", {
  set.seed(2102)
  zi_dat <- tweedie_dat
  zi_dat$y[rbinom(nrow(zi_dat), 1, plogis(-2 + 0.5 * zi_dat$z)) == 1] <- 0

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    ziformula = ~ z,
    family = tweedie(),
    data = zi_dat,
    start = tweedie_start,
    map = tweedie_map,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    ziformula = ~ z,
    family = tweedie(),
    data = zi_dat,
    start = tweedie_start,
    map = tweedie_map,
    se = FALSE
  )

  expect_tweedie_matches(m_rtmb, m_tmb)
  expect_equal(fixef(m_rtmb)$zi, fixef(m_tmb)$zi, tolerance = tol_fixef)
})

test_that("tweedie: prediction and simulation", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    family = tweedie(),
    data = tweedie_dat,
    start = tweedie_start,
    map = tweedie_map,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    family = tweedie(),
    data = tweedie_dat,
    start = tweedie_start,
    map = tweedie_map,
    se = FALSE
  )

  expect_equal(
    predict(m_rtmb, type = "response"),
    predict(m_tmb, type = "response"),
    tolerance = tol_fixef
  )
  expect_equal(
    predict(m_rtmb, type = "disp"),
    predict(m_tmb, type = "disp"),
    tolerance = tol_fixef
  )

  set.seed(2103)
  sim_y <- simulate(m_rtmb, nsim = 1)[[1]]
  expect_equal(length(sim_y), nrow(tweedie_dat))
  expect_true(all(is.finite(sim_y)))
  expect_true(all(sim_y >= 0))
})
