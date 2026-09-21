## Test cases for the RTMB nbinom2 family
## Fit each model with the RTMB and legacy TMB backends, then compare
## likelihoods, fixed effects, and covariance estimates where applicable.

context("RTMB nbinom2 backend")

skip_if_not_installed("RTMB")

data("Salamanders", package = "glmmTMB")

tol_logLik <- 1e-5
tol_fixef <- 1e-5
tol_varcorr <- 1e-4

set.seed(401)
nbinom2_dat <- expand.grid(
  obs = seq_len(5),
  g = factor(seq_len(30))
)
nbinom2_dat$x <- rnorm(nrow(nbinom2_dat))
cond_effect <- rnorm(nlevels(nbinom2_dat$g), sd = 0.35)
zi_effect <- rnorm(nlevels(nbinom2_dat$g), sd = 0.45)
mu <- exp(0.4 + 0.25 * nbinom2_dat$x + cond_effect[nbinom2_dat$g])
zi_prob <- plogis(-1.4 + zi_effect[nbinom2_dat$g])
is_structural_zero <- rbinom(nrow(nbinom2_dat), size = 1, prob = zi_prob)
nbinom2_dat$count <- ifelse(
  is_structural_zero == 1,
  0,
  rnbinom(nrow(nbinom2_dat), mu = mu, size = 2.5)
)
nbinom2_dat$w <- runif(nrow(nbinom2_dat), 0.5, 1.5)
nbinom2_dat$off <- log(runif(nrow(nbinom2_dat), 0.8, 1.2))

set.seed(501)
nbinom2_zi_re_dat <- expand.grid(
  obs = seq_len(8),
  g = factor(seq_len(30))
)
nbinom2_zi_re_dat$x <- rnorm(nrow(nbinom2_zi_re_dat))
zi_re_effect <- rnorm(nlevels(nbinom2_zi_re_dat$g), sd = 0.25)
zi_re_mu <- exp(0.7 + 0.2 * nbinom2_zi_re_dat$x)
zi_re_prob <- plogis(-2.0 + zi_re_effect[nbinom2_zi_re_dat$g])
zi_re_zero <- rbinom(nrow(nbinom2_zi_re_dat), size = 1, prob = zi_re_prob)
nbinom2_zi_re_dat$count <- ifelse(
  zi_re_zero == 1,
  0,
  rnbinom(nrow(nbinom2_zi_re_dat), mu = zi_re_mu, size = 4)
)

test_that("nbinom2: fixed conditional effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ mined,
    family = nbinom2,
    data = Salamanders,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ mined,
    family = nbinom2,
    data = Salamanders,
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
})

test_that("nbinom2: offsets and weights", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x + offset(off),
    weights = w,
    family = nbinom2,
    data = nbinom2_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x + offset(off),
    weights = w,
    family = nbinom2,
    data = nbinom2_dat,
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
})

test_that("nbinom2: dispersion fixed effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x,
    dispformula = ~ x,
    family = nbinom2,
    data = nbinom2_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x,
    dispformula = ~ x,
    family = nbinom2,
    data = nbinom2_dat,
    se = FALSE
  )

  expect_equal(
    as.numeric(logLik(m_rtmb)),
    as.numeric(logLik(m_tmb)),
    tolerance = tol_logLik
  )
  expect_equal(
    fixef(m_rtmb)$disp,
    fixef(m_tmb)$disp,
    tolerance = tol_fixef
  )
})

test_that("nbinom2: conditional random intercept", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x + (1 | g),
    family = nbinom2,
    data = nbinom2_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x + (1 | g),
    family = nbinom2,
    data = nbinom2_dat,
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
    as.numeric(VarCorr(m_rtmb)$cond$g),
    as.numeric(VarCorr(m_tmb)$cond$g),
    tolerance = tol_varcorr
  )
})

test_that("nbinom2: zero-inflation fixed effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x,
    ziformula = ~ x,
    family = nbinom2,
    data = nbinom2_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x,
    ziformula = ~ x,
    family = nbinom2,
    data = nbinom2_dat,
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
})

test_that("nbinom2: zero-inflation random effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x,
    ziformula = ~ 1 + (1 | g),
    family = nbinom2,
    data = nbinom2_zi_re_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x,
    ziformula = ~ 1 + (1 | g),
    family = nbinom2,
    data = nbinom2_zi_re_dat,
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
    as.numeric(VarCorr(m_rtmb)$zi$g),
    as.numeric(VarCorr(m_tmb)$zi$g),
    tolerance = tol_varcorr
  )
})

test_that("nbinom2: prediction with standard errors", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x + (1 | g),
    family = nbinom2,
    data = nbinom2_dat
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x + (1 | g),
    family = nbinom2,
    data = nbinom2_dat
  )

  p_rtmb <- predict(m_rtmb, type = "response", se.fit = TRUE)
  p_tmb <- predict(m_tmb, type = "response", se.fit = TRUE)

  expect_equal(p_rtmb$fit, p_tmb$fit, tolerance = tol_fixef)
  expect_equal(p_rtmb$se.fit, p_tmb$se.fit, tolerance = tol_varcorr)
})

test_that("nbinom2: simulate works under RTMB backend", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x + (1 | g),
    family = nbinom2,
    data = nbinom2_dat,
    se = FALSE
  )
  sim <- m_rtmb$obj$simulate(complete = TRUE)

  expect_true(is.list(sim))
  expect_equal(length(sim$yobs), nrow(nbinom2_dat))
  expect_true(all(sim$yobs >= 0))
  expect_true(all(abs(sim$yobs - round(sim$yobs)) < .Machine$double.eps))
})
