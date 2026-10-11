## Test cases for the RTMB compois family
## Fit each model with the RTMB and legacy TMB backends, then compare
## likelihoods, fixed effects, and covariance estimates where applicable.

context("RTMB compois backend")

skip_if_not_installed("RTMB")

data("Salamanders", package = "glmmTMB")

tol_logLik <- 1e-5
tol_fixef <- 1e-5
tol_varcorr <- 1e-4

set.seed(901)
compois_dat <- expand.grid(
  obs = seq_len(5),
  g = factor(seq_len(28))
)
compois_dat$x <- rnorm(nrow(compois_dat))
cond_effect <- rnorm(nlevels(compois_dat$g), sd = 0.35)
zi_effect <- rnorm(nlevels(compois_dat$g), sd = 0.35)
mu <- exp(0.35 + 0.25 * compois_dat$x + cond_effect[compois_dat$g])
zi_prob <- plogis(-1.6 + zi_effect[compois_dat$g])
is_structural_zero <- rbinom(nrow(compois_dat), size = 1, prob = zi_prob)
compois_dat$count <- ifelse(
  is_structural_zero == 1,
  0,
  rpois(nrow(compois_dat), mu)
)
compois_dat$w <- runif(nrow(compois_dat), 0.5, 1.5)
compois_dat$off <- log(runif(nrow(compois_dat), 0.8, 1.2))

set.seed(902)
compois_zi_re_dat <- expand.grid(
  obs = seq_len(6),
  g = factor(seq_len(25))
)
compois_zi_re_dat$x <- rnorm(nrow(compois_zi_re_dat))
zi_re_effect <- rnorm(nlevels(compois_zi_re_dat$g), sd = 0.3)
zi_re_mu <- exp(0.5 + 0.2 * compois_zi_re_dat$x)
zi_re_prob <- plogis(-1.8 + zi_re_effect[compois_zi_re_dat$g])
zi_re_zero <- rbinom(nrow(compois_zi_re_dat), size = 1, prob = zi_re_prob)
compois_zi_re_dat$count <- ifelse(
  zi_re_zero == 1,
  0,
  rpois(nrow(compois_zi_re_dat), zi_re_mu)
)

test_that("compois: fixed conditional effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ mined,
    family = compois,
    data = Salamanders,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ mined,
    family = compois,
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

test_that("compois: offsets and weights", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x + offset(off),
    weights = w,
    family = compois,
    data = compois_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x + offset(off),
    weights = w,
    family = compois,
    data = compois_dat,
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

test_that("compois: dispersion fixed effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x,
    dispformula = ~ x,
    family = compois,
    data = compois_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x,
    dispformula = ~ x,
    family = compois,
    data = compois_dat,
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

test_that("compois: conditional random intercept", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x + (1 | g),
    family = compois,
    data = compois_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x + (1 | g),
    family = compois,
    data = compois_dat,
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

test_that("compois: zero-inflation fixed effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x,
    ziformula = ~ x,
    family = compois,
    data = compois_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x,
    ziformula = ~ x,
    family = compois,
    data = compois_dat,
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

test_that("compois: zero-inflation random effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x,
    ziformula = ~ 1 + (1 | g),
    family = compois,
    data = compois_zi_re_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x,
    ziformula = ~ 1 + (1 | g),
    family = compois,
    data = compois_zi_re_dat,
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

test_that("compois: prediction with standard errors", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x + (1 | g),
    family = compois,
    data = compois_dat
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x + (1 | g),
    family = compois,
    data = compois_dat
  )

  p_rtmb <- predict(m_rtmb, type = "response", se.fit = TRUE)
  p_tmb <- predict(m_tmb, type = "response", se.fit = TRUE)

  expect_equal(p_rtmb$fit, p_tmb$fit, tolerance = tol_fixef)
  expect_equal(p_rtmb$se.fit, p_tmb$se.fit, tolerance = tol_varcorr)
})

test_that("compois: simulate works under RTMB backend", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x + (1 | g),
    family = compois,
    data = compois_dat,
    se = FALSE
  )
  sim <- m_rtmb$obj$simulate(complete = TRUE)

  expect_true(is.list(sim))
  expect_equal(length(sim$yobs), nrow(compois_dat))
  expect_true(all(sim$yobs >= 0))
  expect_true(all(abs(sim$yobs - round(sim$yobs)) < .Machine$double.eps))
})
