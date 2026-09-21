## Test cases for the RTMB zero-truncated nbinom1 family
## Fit each model with RTMB and legacy TMB, then compare likelihoods,
## fixed effects, and covariance estimates where applicable.

context("RTMB truncated nbinom1 backend")

skip_if_not_installed("RTMB")

data("Salamanders", package = "glmmTMB")

tol_logLik <- 1e-5
tol_fixef <- 1e-5
tol_varcorr <- 1e-4

positive_salamanders <- subset(Salamanders, count > 0)

set.seed(801)
truncated_nbinom1_dat <- expand.grid(
  obs = seq_len(6),
  g = factor(seq_len(25))
)
truncated_nbinom1_dat$x <- rnorm(nrow(truncated_nbinom1_dat))
cond_effect <- rnorm(nlevels(truncated_nbinom1_dat$g), sd = 0.3)
mu <- exp(
  0.7 + 0.2 * truncated_nbinom1_dat$x +
    cond_effect[truncated_nbinom1_dat$g]
)
phi <- 0.8
truncated_nbinom1_dat$count <- glmmTMB:::rtruncated_nbinom_rtmb(
  nrow(truncated_nbinom1_dat),
  size = mu / phi,
  mu = mu
)
truncated_nbinom1_dat$off <- log(
  runif(nrow(truncated_nbinom1_dat), 0.8, 1.2)
)
truncated_nbinom1_dat$w <- runif(nrow(truncated_nbinom1_dat), 0.5, 1.5)

test_that("truncated nbinom1: RTMB density matches package density", {
  x <- 1:8
  mu <- seq(0.5, 4, length.out = length(x))
  phi <- seq(0.4, 1.2, length.out = length(x))
  log_mu <- log(mu)
  log_phi <- log(phi)
  log_var_minus_mu <- log_mu + log_phi

  expect_equal(
    glmmTMB:::dtruncated_nbinom1_rtmb(
      x, log_mu, log_var_minus_mu, log_phi, log = TRUE
    ),
    glmmTMB::dtruncated_nbinom1(x, phi = phi, mu = mu, log = TRUE),
    tolerance = 1e-12
  )
  expect_equal(
    glmmTMB:::dtruncated_nbinom1_rtmb(
      0, log(1), log(1), log(1), log = TRUE
    ),
    -Inf
  )
})

test_that("truncated nbinom1: fixed conditional effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ mined,
    family = truncated_nbinom1,
    data = positive_salamanders,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ mined,
    family = truncated_nbinom1,
    data = positive_salamanders,
    se = FALSE
  )

  expect_equal(
    as.numeric(logLik(m_rtmb)),
    as.numeric(logLik(m_tmb)),
    tolerance = tol_logLik
  )
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = tol_fixef)
  expect_equal(sigma(m_rtmb), sigma(m_tmb), tolerance = tol_fixef)
})

test_that("truncated nbinom1: offsets and weights", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x + offset(off),
    weights = w,
    family = truncated_nbinom1,
    data = truncated_nbinom1_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x + offset(off),
    weights = w,
    family = truncated_nbinom1,
    data = truncated_nbinom1_dat,
    se = FALSE
  )

  expect_equal(
    as.numeric(logLik(m_rtmb)),
    as.numeric(logLik(m_tmb)),
    tolerance = tol_logLik
  )
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = tol_fixef)
})

test_that("truncated nbinom1: dispersion fixed effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x,
    dispformula = ~ x,
    family = truncated_nbinom1,
    data = truncated_nbinom1_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x,
    dispformula = ~ x,
    family = truncated_nbinom1,
    data = truncated_nbinom1_dat,
    se = FALSE
  )

  expect_equal(
    as.numeric(logLik(m_rtmb)),
    as.numeric(logLik(m_tmb)),
    tolerance = tol_logLik
  )
  expect_equal(fixef(m_rtmb)$disp, fixef(m_tmb)$disp, tolerance = tol_fixef)
})

test_that("truncated nbinom1: conditional random intercept", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x + (1 | g),
    family = truncated_nbinom1,
    data = truncated_nbinom1_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x + (1 | g),
    family = truncated_nbinom1,
    data = truncated_nbinom1_dat,
    se = FALSE
  )

  expect_equal(
    as.numeric(logLik(m_rtmb)),
    as.numeric(logLik(m_tmb)),
    tolerance = tol_logLik
  )
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = tol_fixef)
  expect_equal(
    as.numeric(VarCorr(m_rtmb)$cond$g),
    as.numeric(VarCorr(m_tmb)$cond$g),
    tolerance = tol_varcorr
  )
})

test_that("truncated nbinom1: fixed zero inflation", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ mined,
    ziformula = ~ mined,
    family = truncated_nbinom1,
    data = Salamanders,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ mined,
    ziformula = ~ mined,
    family = truncated_nbinom1,
    data = Salamanders,
    se = FALSE
  )

  expect_equal(
    as.numeric(logLik(m_rtmb)),
    as.numeric(logLik(m_tmb)),
    tolerance = tol_logLik
  )
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = tol_fixef)
  expect_equal(fixef(m_rtmb)$zi, fixef(m_tmb)$zi, tolerance = tol_fixef)
})

test_that("truncated nbinom1: prediction with standard errors", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x + (1 | g),
    family = truncated_nbinom1,
    data = truncated_nbinom1_dat
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x + (1 | g),
    family = truncated_nbinom1,
    data = truncated_nbinom1_dat
  )

  p_rtmb <- predict(m_rtmb, type = "response", se.fit = TRUE)
  p_tmb <- predict(m_tmb, type = "response", se.fit = TRUE)

  expect_equal(p_rtmb$fit, p_tmb$fit, tolerance = tol_fixef)
  expect_equal(p_rtmb$se.fit, p_tmb$se.fit, tolerance = tol_varcorr)
})

test_that("truncated nbinom1: simulation is strictly positive", {
  local_useRTMB(TRUE)
  model <- glmmTMB(
    count ~ x + (1 | g),
    family = truncated_nbinom1,
    data = truncated_nbinom1_dat,
    se = FALSE
  )

  set.seed(802)
  simulated <- model$obj$simulate(complete = TRUE)$yobs

  expect_length(simulated, nrow(truncated_nbinom1_dat))
  expect_true(all(is.finite(simulated)))
  expect_true(all(simulated > 0))
  expect_true(all(simulated == floor(simulated)))
})

