## Test cases for the RTMB binomial family
## Fit each model with the RTMB and legacy TMB backends, then compare
## likelihoods, fixed effects, and covariance estimates where applicable.

context("RTMB Binomial backend")

skip_if_not_installed("RTMB")

data("cbpp", package = "lme4")

tol_logLik <- 1e-5
tol_fixef <- 1e-5
tol_varcorr <- 1e-4

set.seed(301)
binom_dat <- data.frame(
  y = rbinom(120, size = 1, prob = 0.35),
  x = rnorm(120),
  g = factor(rep(seq_len(24), each = 5))
)

set.seed(302)
zi_binom_dat <- expand.grid(
  obs = seq_len(4),
  g = factor(seq_len(20))
)
zi_binom_dat$x <- rnorm(nrow(zi_binom_dat))
zi_binom_dat$size <- 10
zi_effect <- rnorm(nlevels(zi_binom_dat$g), sd = 0.5)
zi_prob <- plogis(-1.2 + zi_effect[zi_binom_dat$g])
cond_prob <- plogis(-0.4 + 0.3 * zi_binom_dat$x)
is_structural_zero <- rbinom(nrow(zi_binom_dat), size = 1, prob = zi_prob)
zi_binom_dat$y <- ifelse(
  is_structural_zero == 1,
  0,
  rbinom(nrow(zi_binom_dat), size = zi_binom_dat$size, prob = cond_prob)
)

test_that("binomial: binary response with fixed effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    family = binomial,
    data = binom_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    family = binomial,
    data = binom_dat,
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

test_that("binomial: empty fixed-effect prior range is skipped", {
  local_useRTMB(TRUE)
  m0 <- glmmTMB(
    y ~ 1,
    family = binomial,
    data = binom_dat,
    se = FALSE
  )
  m_prior <- expect_no_error(glmmTMB(
    y ~ 1,
    family = binomial,
    data = binom_dat,
    priors = data.frame(prior = "normal(0, 3)", class = "beta"),
    se = FALSE
  ))

  expect_equal(logLik(m_prior), logLik(m0), tolerance = tol_logLik)
})

test_that("binomial: grouped response with weights", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    incidence / size ~ period,
    weights = size,
    family = binomial,
    data = cbpp,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    incidence / size ~ period,
    weights = size,
    family = binomial,
    data = cbpp,
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

test_that("binomial: cbind response", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    cbind(incidence, size - incidence) ~ period,
    family = binomial,
    data = cbpp,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    cbind(incidence, size - incidence) ~ period,
    family = binomial,
    data = cbpp,
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

test_that("binomial: conditional random intercept", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    incidence / size ~ period + (1 | herd),
    weights = size,
    family = binomial,
    data = cbpp,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    incidence / size ~ period + (1 | herd),
    weights = size,
    family = binomial,
    data = cbpp,
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
    as.numeric(VarCorr(m_rtmb)$cond$herd),
    as.numeric(VarCorr(m_tmb)$cond$herd),
    tolerance = tol_varcorr
  )
})

test_that("binomial: cloglog link", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    family = binomial(link = "cloglog"),
    data = binom_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    family = binomial(link = "cloglog"),
    data = binom_dat,
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

test_that("binomial: zero-inflation fixed effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    ziformula = ~ x,
    family = binomial,
    data = binom_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    ziformula = ~ x,
    family = binomial,
    data = binom_dat,
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
    fixef(m_rtmb)$zi,
    fixef(m_tmb)$zi,
    tolerance = tol_fixef
  )
})

test_that("binomial: zero-inflation random effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y / size ~ x,
    ziformula = ~ 1 + (1 | g),
    weights = size,
    family = binomial,
    data = zi_binom_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y / size ~ x,
    ziformula = ~ 1 + (1 | g),
    weights = size,
    family = binomial,
    data = zi_binom_dat,
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

test_that("binomial: simulate works under RTMB backend", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    incidence / size ~ period + (1 | herd),
    weights = size,
    family = binomial,
    data = cbpp,
    se = FALSE
  )
  sim <- m_rtmb$obj$simulate(complete = TRUE)

  expect_true(is.list(sim))
  expect_equal(length(sim$yobs), nrow(cbpp))
  expect_true(all(sim$yobs >= 0))
  expect_true(all(sim$yobs <= cbpp$size))
})
