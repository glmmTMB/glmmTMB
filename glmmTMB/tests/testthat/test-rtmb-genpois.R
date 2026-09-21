## Test cases for the RTMB generalized Poisson family
## Fit each model with the RTMB and legacy TMB backends, then compare
## likelihoods, fixed effects, and covariance estimates where applicable.

context("RTMB genpois backend")

skip_if_not_installed("RTMB")

tol_logLik <- 1e-5
tol_fixef <- 1e-5
tol_varcorr <- 1e-4

set.seed(1401)
genpois_dat <- expand.grid(
  obs = seq_len(5),
  g = factor(seq_len(30))
)
genpois_dat$x <- rnorm(nrow(genpois_dat))
cond_effect <- rnorm(nlevels(genpois_dat$g), sd = 0.25)
zi_effect <- rnorm(nlevels(genpois_dat$g), sd = 0.3)
mu <- exp(0.6 + 0.25 * genpois_dat$x + cond_effect[genpois_dat$g])
phi <- 1.5
theta <- mu / sqrt(phi)
lambda <- 1 - 1 / sqrt(phi)
zi_prob <- plogis(-1.7 + zi_effect[genpois_dat$g])
is_structural_zero <- rbinom(nrow(genpois_dat), size = 1, prob = zi_prob)
genpois_dat$count <- ifelse(
  is_structural_zero == 1,
  0,
  glmmTMB::rgenpois(nrow(genpois_dat), theta, lambda)
)
genpois_dat$w <- runif(nrow(genpois_dat), 0.5, 1.5)
genpois_dat$off <- log(runif(nrow(genpois_dat), 0.8, 1.2))

expect_genpois_equal <- function(m_rtmb, m_tmb, check_varcorr = FALSE) {
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
    fixef(m_rtmb)$disp,
    fixef(m_tmb)$disp,
    tolerance = tol_fixef
  )
  expect_equal(
    fixef(m_rtmb)$zi,
    fixef(m_tmb)$zi,
    tolerance = tol_fixef
  )
  if (check_varcorr) {
    expect_equal(
      unname(unlist(VarCorr(m_rtmb)$cond)),
      unname(unlist(VarCorr(m_tmb)$cond)),
      tolerance = tol_varcorr
    )
  }
}

test_that("genpois: RTMB density matches package density", {
  x <- 0:8
  theta <- 1.4
  lambda <- 0.2

  expect_equal(
    glmmTMB:::dgenpois_rtmb(x, theta = theta, lambda = lambda, log = TRUE),
    glmmTMB::dgenpois(x, lambda1 = theta, lambda2 = lambda, log = TRUE),
    tolerance = 1e-12
  )
  expect_equal(
    glmmTMB:::dgenpois_rtmb(x, theta = theta, lambda = lambda, log = FALSE),
    glmmTMB::dgenpois(x, lambda1 = theta, lambda2 = lambda, log = FALSE),
    tolerance = 1e-12
  )
})

test_that("genpois: fixed conditional effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x,
    family = genpois,
    data = genpois_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x,
    family = genpois,
    data = genpois_dat,
    se = FALSE
  )

  expect_genpois_equal(m_rtmb, m_tmb)
})

test_that("genpois: offsets and weights", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x + offset(off),
    weights = w,
    family = genpois,
    data = genpois_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x + offset(off),
    weights = w,
    family = genpois,
    data = genpois_dat,
    se = FALSE
  )

  expect_genpois_equal(m_rtmb, m_tmb)
})

test_that("genpois: dispersion fixed effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x,
    dispformula = ~ x,
    family = genpois,
    data = genpois_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x,
    dispformula = ~ x,
    family = genpois,
    data = genpois_dat,
    se = FALSE
  )

  expect_genpois_equal(m_rtmb, m_tmb)
})

test_that("genpois: conditional random intercept", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x + (1 | g),
    family = genpois,
    data = genpois_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x + (1 | g),
    family = genpois,
    data = genpois_dat,
    se = FALSE
  )

  expect_genpois_equal(m_rtmb, m_tmb, check_varcorr = TRUE)
})

test_that("genpois: zero-inflation fixed effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x,
    ziformula = ~ x,
    family = genpois,
    data = genpois_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x,
    ziformula = ~ x,
    family = genpois,
    data = genpois_dat,
    se = FALSE
  )

  expect_genpois_equal(m_rtmb, m_tmb)
})

test_that("genpois: prediction with standard errors", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x,
    family = genpois,
    data = genpois_dat
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x,
    family = genpois,
    data = genpois_dat
  )

  expect_equal(
    predict(m_rtmb, type = "response"),
    predict(m_tmb, type = "response"),
    tolerance = tol_fixef
  )
  expect_equal(
    predict(m_rtmb, type = "response", se.fit = TRUE)$fit,
    predict(m_tmb, type = "response", se.fit = TRUE)$fit,
    tolerance = tol_fixef
  )
})

test_that("genpois: simulate works under RTMB backend", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x,
    family = genpois,
    data = genpois_dat,
    se = FALSE
  )

  set.seed(1402)
  sim <- simulate(m_rtmb, nsim = 1)
  simulated <- sim[[1]]

  expect_length(simulated, nrow(genpois_dat))
  expect_true(all(simulated >= 0))
  expect_true(all(simulated == floor(simulated)))
})
