## Test cases for the RTMB nbinom1 family
## Fit each model with the RTMB and legacy TMB backends, then compare
## likelihoods, fixed effects, and covariance estimates where applicable.

context("RTMB nbinom1 backend")

skip_if_not_installed("RTMB")

data("Salamanders", package = "glmmTMB")

tol_logLik <- 1e-5
tol_fixef <- 1e-5
tol_varcorr <- 1e-4

set.seed(411)
nbinom1_dat <- expand.grid(
  obs = seq_len(5),
  g = factor(seq_len(30))
)
nbinom1_dat$x <- rnorm(nrow(nbinom1_dat))
cond_effect <- rnorm(nlevels(nbinom1_dat$g), sd = 0.35)
zi_effect <- rnorm(nlevels(nbinom1_dat$g), sd = 0.45)
mu <- exp(0.4 + 0.25 * nbinom1_dat$x + cond_effect[nbinom1_dat$g])
phi <- 1.4
zi_prob <- plogis(-1.4 + zi_effect[nbinom1_dat$g])
is_structural_zero <- rbinom(nrow(nbinom1_dat), size = 1, prob = zi_prob)
nbinom1_dat$count <- ifelse(
  is_structural_zero == 1,
  0,
  rnbinom(nrow(nbinom1_dat), mu = mu, size = mu / phi)
)
nbinom1_dat$w <- runif(nrow(nbinom1_dat), 0.5, 1.5)
nbinom1_dat$off <- log(runif(nrow(nbinom1_dat), 0.8, 1.2))

set.seed(601)
nbinom1_zi_re_dat <- expand.grid(
  obs = seq_len(10),
  g = factor(seq_len(35))
)
nbinom1_zi_re_dat$x <- rnorm(nrow(nbinom1_zi_re_dat))
zi_re_effect <- rnorm(nlevels(nbinom1_zi_re_dat$g), sd = 0.2)
zi_re_mu <- exp(0.8 + 0.15 * nbinom1_zi_re_dat$x)
zi_re_phi <- 0.8
zi_re_prob <- plogis(-2.2 + zi_re_effect[nbinom1_zi_re_dat$g])
zi_re_zero <- rbinom(nrow(nbinom1_zi_re_dat), size = 1, prob = zi_re_prob)
nbinom1_zi_re_dat$count <- ifelse(
  zi_re_zero == 1,
  0,
  rnbinom(nrow(nbinom1_zi_re_dat), mu = zi_re_mu,
          size = zi_re_mu / zi_re_phi)
)

test_that("nbinom1: fixed conditional effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ mined,
    family = nbinom1,
    data = Salamanders,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ mined,
    family = nbinom1,
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

test_that("nbinom1: offsets and weights", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x + offset(off),
    weights = w,
    family = nbinom1,
    data = nbinom1_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x + offset(off),
    weights = w,
    family = nbinom1,
    data = nbinom1_dat,
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

test_that("nbinom1: dispersion fixed effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x,
    dispformula = ~ x,
    family = nbinom1,
    data = nbinom1_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x,
    dispformula = ~ x,
    family = nbinom1,
    data = nbinom1_dat,
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

test_that("nbinom1: conditional random intercept", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x + (1 | g),
    family = nbinom1,
    data = nbinom1_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x + (1 | g),
    family = nbinom1,
    data = nbinom1_dat,
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

test_that("nbinom1: zero-inflation fixed effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x,
    ziformula = ~ x,
    family = nbinom1,
    data = nbinom1_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x,
    ziformula = ~ x,
    family = nbinom1,
    data = nbinom1_dat,
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

test_that("nbinom1: zero-inflation random effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x,
    ziformula = ~ 1 + (1 | g),
    family = nbinom1,
    data = nbinom1_zi_re_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x,
    ziformula = ~ 1 + (1 | g),
    family = nbinom1,
    data = nbinom1_zi_re_dat,
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

test_that("nbinom1: prediction with standard errors", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x + (1 | g),
    family = nbinom1,
    data = nbinom1_dat
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x + (1 | g),
    family = nbinom1,
    data = nbinom1_dat
  )

  p_rtmb <- predict(m_rtmb, type = "response", se.fit = TRUE)
  p_tmb <- predict(m_tmb, type = "response", se.fit = TRUE)

  expect_equal(p_rtmb$fit, p_tmb$fit, tolerance = tol_fixef)
  expect_equal(p_rtmb$se.fit, p_tmb$se.fit, tolerance = tol_varcorr)
})

test_that("nbinom1: simulate works under RTMB backend", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x + (1 | g),
    family = nbinom1,
    data = nbinom1_dat,
    se = FALSE
  )
  sim <- m_rtmb$obj$simulate(complete = TRUE)

  expect_true(is.list(sim))
  expect_equal(length(sim$yobs), nrow(nbinom1_dat))
  expect_true(all(sim$yobs >= 0))
  expect_true(all(abs(sim$yobs - round(sim$yobs)) < .Machine$double.eps))
})
