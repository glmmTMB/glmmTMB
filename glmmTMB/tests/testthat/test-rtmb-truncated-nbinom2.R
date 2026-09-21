## Test cases for the RTMB zero-truncated nbinom2 family
## Fit each model with RTMB and legacy TMB, then compare likelihoods,
## fixed effects, and covariance estimates where applicable.

context("RTMB truncated nbinom2 backend")

skip_if_not_installed("RTMB")

data("Salamanders", package = "glmmTMB")

tol_logLik <- 1e-5
tol_fixef <- 1e-5
tol_varcorr <- 1e-4

positive_salamanders <- subset(Salamanders, count > 0)

set.seed(701)
truncated_nbinom2_dat <- expand.grid(
  obs = seq_len(6),
  g = factor(seq_len(25))
)
truncated_nbinom2_dat$x <- rnorm(nrow(truncated_nbinom2_dat))
cond_effect <- rnorm(nlevels(truncated_nbinom2_dat$g), sd = 0.3)
mu <- exp(
  0.7 + 0.2 * truncated_nbinom2_dat$x +
    cond_effect[truncated_nbinom2_dat$g]
)
truncated_nbinom2_dat$count <- glmmTMB:::rtruncated_nbinom_rtmb(
  nrow(truncated_nbinom2_dat),
  size = 3,
  mu = mu
)
truncated_nbinom2_dat$off <- log(
  runif(nrow(truncated_nbinom2_dat), 0.8, 1.2)
)
truncated_nbinom2_dat$w <- runif(nrow(truncated_nbinom2_dat), 0.5, 1.5)

test_that("truncated nbinom2: RTMB density matches package density", {
  x <- 1:8
  mu <- seq(0.5, 4, length.out = length(x))
  size <- seq(1, 3, length.out = length(x))
  log_mu <- log(mu)
  log_size <- log(size)
  log_var_minus_mu <- 2 * log_mu - log_size

  expect_equal(
    glmmTMB:::dtruncated_nbinom2_rtmb(
      x, log_mu, log_var_minus_mu, log_size, log = TRUE
    ),
    glmmTMB::dtruncated_nbinom2(x, size = size, mu = mu, log = TRUE),
    tolerance = 1e-12
  )
  expect_equal(
    glmmTMB:::dtruncated_nbinom2_rtmb(
      0, log(1), log(1), log(1), log = TRUE
    ),
    -Inf
  )
})

test_that("truncated nbinom2: fixed conditional effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ mined,
    family = truncated_nbinom2,
    data = positive_salamanders,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ mined,
    family = truncated_nbinom2,
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

test_that("truncated nbinom2: offsets and weights", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x + offset(off),
    weights = w,
    family = truncated_nbinom2,
    data = truncated_nbinom2_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x + offset(off),
    weights = w,
    family = truncated_nbinom2,
    data = truncated_nbinom2_dat,
    se = FALSE
  )

  expect_equal(
    as.numeric(logLik(m_rtmb)),
    as.numeric(logLik(m_tmb)),
    tolerance = tol_logLik
  )
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = tol_fixef)
})

test_that("truncated nbinom2: dispersion fixed effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x,
    dispformula = ~ x,
    family = truncated_nbinom2,
    data = truncated_nbinom2_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x,
    dispformula = ~ x,
    family = truncated_nbinom2,
    data = truncated_nbinom2_dat,
    se = FALSE
  )

  expect_equal(
    as.numeric(logLik(m_rtmb)),
    as.numeric(logLik(m_tmb)),
    tolerance = tol_logLik
  )
  expect_equal(fixef(m_rtmb)$disp, fixef(m_tmb)$disp, tolerance = tol_fixef)
})

test_that("truncated nbinom2: conditional random intercept", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x + (1 | g),
    family = truncated_nbinom2,
    data = truncated_nbinom2_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x + (1 | g),
    family = truncated_nbinom2,
    data = truncated_nbinom2_dat,
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

test_that("truncated nbinom2: fixed zero inflation", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ mined,
    ziformula = ~ mined,
    family = truncated_nbinom2,
    data = Salamanders,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ mined,
    ziformula = ~ mined,
    family = truncated_nbinom2,
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

test_that("truncated nbinom2: zero-inflation random intercept", {
  thetazi <- log(0.5)
  theta_map <- factor(NA)

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ mined,
    ziformula = ~ mined + (1 | site),
    family = truncated_nbinom2,
    data = Salamanders,
    start = list(thetazi = thetazi),
    map = list(thetazi = theta_map),
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ mined,
    ziformula = ~ mined + (1 | site),
    family = truncated_nbinom2,
    data = Salamanders,
    start = list(thetazi = thetazi),
    map = list(thetazi = theta_map),
    se = FALSE
  )

  expect_equal(
    as.numeric(logLik(m_rtmb)),
    as.numeric(logLik(m_tmb)),
    tolerance = tol_logLik
  )
  expect_equal(fixef(m_rtmb)$zi, fixef(m_tmb)$zi, tolerance = tol_fixef)
  expect_equal(
    as.numeric(VarCorr(m_rtmb)$zi$site),
    as.numeric(VarCorr(m_tmb)$zi$site),
    tolerance = tol_varcorr
  )
})

test_that("truncated nbinom2: prediction with standard errors", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x + (1 | g),
    family = truncated_nbinom2,
    data = truncated_nbinom2_dat
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x + (1 | g),
    family = truncated_nbinom2,
    data = truncated_nbinom2_dat
  )

  p_rtmb <- predict(m_rtmb, type = "response", se.fit = TRUE)
  p_tmb <- predict(m_tmb, type = "response", se.fit = TRUE)

  expect_equal(p_rtmb$fit, p_tmb$fit, tolerance = tol_fixef)
  expect_equal(p_rtmb$se.fit, p_tmb$se.fit, tolerance = tol_varcorr)
})

test_that("truncated nbinom2: simulation is strictly positive", {
  local_useRTMB(TRUE)
  model <- glmmTMB(
    count ~ x + (1 | g),
    family = truncated_nbinom2,
    data = truncated_nbinom2_dat,
    se = FALSE
  )

  set.seed(702)
  simulated <- model$obj$simulate(complete = TRUE)$yobs

  expect_length(simulated, nrow(truncated_nbinom2_dat))
  expect_true(all(is.finite(simulated)))
  expect_true(all(simulated > 0))
  expect_true(all(simulated == floor(simulated)))
})

test_that("truncated nbinom2: ZI simulation contains zeros and positives", {
  fixed_map <- factor(NA)

  local_useRTMB(TRUE)
  model <- glmmTMB(
    count ~ 1,
    ziformula = ~ 1,
    family = truncated_nbinom2,
    data = Salamanders,
    start = list(
      beta = log(5),
      betadisp = log(2),
      betazi = qlogis(0.5)
    ),
    map = list(
      beta = fixed_map,
      betadisp = fixed_map,
      betazi = fixed_map
    ),
    se = FALSE
  )

  set.seed(703)
  simulated <- model$obj$simulate(complete = TRUE)$yobs

  expect_length(simulated, nrow(Salamanders))
  expect_true(all(is.finite(simulated)))
  expect_true(any(simulated == 0))
  expect_true(any(simulated > 0))
  expect_true(all(simulated == floor(simulated)))
})
