## Test cases for the RTMB nbinom12 family
## Fit each model with RTMB and legacy TMB, then compare likelihoods,
## fixed effects, dispersion, and the extra nbinom12 psi parameter.

context("RTMB nbinom12 backend")

skip_if_not_installed("RTMB")

data("Salamanders", package = "glmmTMB")

tol_logLik <- 1e-5
tol_fixef <- 1e-5
tol_varcorr <- 1e-4

set.seed(901)
nbinom12_dat <- expand.grid(
  obs = seq_len(6),
  g = factor(seq_len(25))
)
nbinom12_dat$x <- rnorm(nrow(nbinom12_dat))
cond_effect <- rnorm(nlevels(nbinom12_dat$g), sd = 0.25)
mu <- exp(
  0.6 + 0.25 * nbinom12_dat$x +
    cond_effect[nbinom12_dat$g]
)
phi <- 0.8
psi <- 2
nbinom12_dat$count <- rnbinom(
  nrow(nbinom12_dat),
  mu = mu,
  size = mu / (phi + mu / psi)
)
nbinom12_dat$off <- log(runif(nrow(nbinom12_dat), 0.8, 1.2))
nbinom12_dat$w <- runif(nrow(nbinom12_dat), 0.5, 1.5)

expect_nbinom12_equal <- function(m_rtmb, m_tmb) {
  expect_equal(
    as.numeric(logLik(m_rtmb)),
    as.numeric(logLik(m_tmb)),
    tolerance = tol_logLik
  )
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = tol_fixef)
  expect_equal(
    m_rtmb$obj$env$last.par.best["psi"],
    m_tmb$obj$env$last.par.best["psi"],
    tolerance = tol_fixef
  )
}

test_that("nbinom12: robust density formula matches stats::dnbinom", {
  x <- 0:8
  mu <- seq(0.5, 4, length.out = length(x))
  phi <- seq(0.4, 1.2, length.out = length(x))
  psi <- 2
  log_mu <- log(mu)
  log_var_minus_mu <- log_mu + log(phi + mu / psi)
  size <- mu / (phi + mu / psi)

  expect_equal(
    glmmTMB:::dnbinom_robust_rtmb(
      x, log_mu = log_mu, log_var_minus_mu = log_var_minus_mu, log = TRUE
    ),
    dnbinom(x, mu = mu, size = size, log = TRUE),
    tolerance = 1e-12
  )
})

test_that("nbinom12: fixed conditional effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x,
    family = nbinom12,
    data = nbinom12_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x,
    family = nbinom12,
    data = nbinom12_dat,
    se = FALSE
  )

  expect_nbinom12_equal(m_rtmb, m_tmb)
  expect_equal(sigma(m_rtmb), sigma(m_tmb), tolerance = tol_fixef)
})

test_that("nbinom12: offsets and weights", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x + offset(off),
    weights = w,
    family = nbinom12,
    data = nbinom12_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x + offset(off),
    weights = w,
    family = nbinom12,
    data = nbinom12_dat,
    se = FALSE
  )

  expect_nbinom12_equal(m_rtmb, m_tmb)
})

test_that("nbinom12: dispersion fixed effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x,
    dispformula = ~ x,
    family = nbinom12,
    data = nbinom12_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x,
    dispformula = ~ x,
    family = nbinom12,
    data = nbinom12_dat,
    se = FALSE
  )

  expect_nbinom12_equal(m_rtmb, m_tmb)
  expect_equal(fixef(m_rtmb)$disp, fixef(m_tmb)$disp, tolerance = tol_fixef)
})

test_that("nbinom12: conditional random intercept", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x + (1 | g),
    family = nbinom12,
    data = nbinom12_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x + (1 | g),
    family = nbinom12,
    data = nbinom12_dat,
    se = FALSE
  )

  expect_nbinom12_equal(m_rtmb, m_tmb)
  expect_equal(
    as.numeric(VarCorr(m_rtmb)$cond$g),
    as.numeric(VarCorr(m_tmb)$cond$g),
    tolerance = tol_varcorr
  )
})

test_that("nbinom12: fixed zero inflation", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ mined,
    ziformula = ~ mined,
    family = nbinom12,
    data = Salamanders,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ mined,
    ziformula = ~ mined,
    family = nbinom12,
    data = Salamanders,
    se = FALSE
  )

  expect_nbinom12_equal(m_rtmb, m_tmb)
  expect_equal(fixef(m_rtmb)$zi, fixef(m_tmb)$zi, tolerance = tol_fixef)
})

test_that("nbinom12: prediction with standard errors", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x + (1 | g),
    family = nbinom12,
    data = nbinom12_dat
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x + (1 | g),
    family = nbinom12,
    data = nbinom12_dat
  )

  p_rtmb <- predict(m_rtmb, type = "response", se.fit = TRUE)
  p_tmb <- predict(m_tmb, type = "response", se.fit = TRUE)

  expect_equal(p_rtmb$fit, p_tmb$fit, tolerance = tol_fixef)
  expect_equal(p_rtmb$se.fit, p_tmb$se.fit, tolerance = tol_varcorr)
})

test_that("nbinom12: simulation works under RTMB backend", {
  local_useRTMB(TRUE)
  model <- glmmTMB(
    count ~ x + (1 | g),
    family = nbinom12,
    data = nbinom12_dat,
    se = FALSE
  )

  set.seed(902)
  simulated <- model$obj$simulate(complete = TRUE)$yobs

  expect_length(simulated, nrow(nbinom12_dat))
  expect_true(all(is.finite(simulated)))
  expect_true(all(simulated >= 0))
  expect_true(all(simulated == floor(simulated)))
})

