## Test cases for the RTMB ordered beta family
## Fit each model with RTMB and legacy TMB, then compare likelihoods,
## parameters, predictions, and covariance estimates where applicable.

context("RTMB ordered beta backend")

skip_if_not_installed("RTMB")

tol_logLik <- 1e-5
tol_fixef <- 1e-5
tol_varcorr <- 1e-4

set.seed(2001)
ordbeta_dat <- expand.grid(
  obs = seq_len(5),
  g = factor(seq_len(30))
)
ordbeta_dat$x <- rnorm(nrow(ordbeta_dat))
ordbeta_dat$z <- rnorm(nrow(ordbeta_dat))
ordbeta_dat$off <- runif(nrow(ordbeta_dat), -0.15, 0.15)
ordbeta_dat$w <- runif(nrow(ordbeta_dat), 0.5, 1.5)
ordbeta_effect <- rnorm(nlevels(ordbeta_dat$g), sd = 0.25)
ordbeta_eta <- -0.2 + 0.35 * ordbeta_dat$x + ordbeta_effect[ordbeta_dat$g]
ordbeta_mu <- plogis(ordbeta_eta)
ordbeta_phi <- 8
ordbeta_cutpoints <- c(-1, 1)
ordbeta_u <- runif(nrow(ordbeta_dat))
ordbeta_p0 <- plogis(ordbeta_cutpoints[1L] - ordbeta_eta)
ordbeta_p1 <- plogis(ordbeta_eta - ordbeta_cutpoints[2L])
ordbeta_dat$y <- ifelse(
  ordbeta_u < ordbeta_p0,
  0,
  ifelse(
    ordbeta_u > 1 - ordbeta_p1,
    1,
    rbeta(
      nrow(ordbeta_dat),
      shape1 = ordbeta_mu * ordbeta_phi,
      shape2 = (1 - ordbeta_mu) * ordbeta_phi
    )
  )
)

ordbeta_start <- list(psi = ordbeta_cutpoints)
ordbeta_map <- list(psi = factor(c(NA, NA)))

expect_ordbeta_matches <- function(m_rtmb, m_tmb, check_varcorr = FALSE) {
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

test_that("ordbeta: RTMB density matches C++ parameterization", {
  x <- c(0, 0.25, 0.75, 1)
  eta <- c(-0.5, -0.2, 0.2, 0.5)
  mu <- plogis(eta)
  phi <- rep(6, length(x))
  cutpoints <- c(-1, 1)

  expected <- c(
    log(plogis(cutpoints[1L] - eta[1L])),
    log(plogis(eta[2L] - cutpoints[1L]) -
          plogis(eta[2L] - cutpoints[2L])) +
      dbeta(x[2L], mu[2L] * phi[2L], (1 - mu[2L]) * phi[2L], log = TRUE),
    log(plogis(eta[3L] - cutpoints[1L]) -
          plogis(eta[3L] - cutpoints[2L])) +
      dbeta(x[3L], mu[3L] * phi[3L], (1 - mu[3L]) * phi[3L], log = TRUE),
    log(plogis(eta[4L] - cutpoints[2L]))
  )

  expect_equal(
    glmmTMB:::dordbeta_rtmb(
      x, eta = eta, mean = mu, phi = phi, cutpoints = cutpoints,
      log = TRUE
    ),
    expected,
    tolerance = 1e-12
  )
})

test_that("ordbeta: fixed conditional effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    family = ordbeta(),
    data = ordbeta_dat,
    start = ordbeta_start,
    map = ordbeta_map,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    family = ordbeta(),
    data = ordbeta_dat,
    start = ordbeta_start,
    map = ordbeta_map,
    se = FALSE
  )

  expect_ordbeta_matches(m_rtmb, m_tmb)
})

test_that("ordbeta: offsets, weights, and dispersion effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x + offset(off),
    dispformula = ~ z,
    weights = w,
    family = ordbeta(),
    data = ordbeta_dat,
    start = ordbeta_start,
    map = ordbeta_map,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x + offset(off),
    dispformula = ~ z,
    weights = w,
    family = ordbeta(),
    data = ordbeta_dat,
    start = ordbeta_start,
    map = ordbeta_map,
    se = FALSE
  )

  expect_ordbeta_matches(m_rtmb, m_tmb)
  expect_equal(fixef(m_rtmb)$disp, fixef(m_tmb)$disp, tolerance = tol_fixef)
})

test_that("ordbeta: conditional random intercept", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x + (1 | g),
    family = ordbeta(),
    data = ordbeta_dat,
    start = ordbeta_start,
    map = ordbeta_map,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x + (1 | g),
    family = ordbeta(),
    data = ordbeta_dat,
    start = ordbeta_start,
    map = ordbeta_map,
    se = FALSE
  )

  expect_ordbeta_matches(m_rtmb, m_tmb, check_varcorr = TRUE)
})

test_that("ordbeta: estimated cutpoints", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    family = ordbeta(),
    data = ordbeta_dat,
    start = ordbeta_start,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    family = ordbeta(),
    data = ordbeta_dat,
    start = ordbeta_start,
    se = FALSE
  )

  expect_ordbeta_matches(m_rtmb, m_tmb)
  rtmb_psi <- m_rtmb$obj$env$last.par.best[
    grepl("^psi", names(m_rtmb$obj$env$last.par.best))
  ]
  tmb_psi <- m_tmb$obj$env$last.par.best[
    grepl("^psi", names(m_tmb$obj$env$last.par.best))
  ]
  expect_equal(
    unname(rtmb_psi),
    unname(tmb_psi),
    tolerance = tol_fixef
  )
  expect_equal(
    family_params(m_rtmb),
    family_params(m_tmb),
    tolerance = tol_fixef
  )
})

test_that("ordbeta: zero inflation", {
  set.seed(2002)
  zi_dat <- ordbeta_dat
  zi_dat$y[rbinom(nrow(zi_dat), 1, plogis(-2 + 0.5 * zi_dat$z)) == 1] <- 0

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    ziformula = ~ z,
    family = ordbeta(),
    data = zi_dat,
    start = ordbeta_start,
    map = ordbeta_map,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    ziformula = ~ z,
    family = ordbeta(),
    data = zi_dat,
    start = ordbeta_start,
    map = ordbeta_map,
    se = FALSE
  )

  expect_ordbeta_matches(m_rtmb, m_tmb)
  expect_equal(fixef(m_rtmb)$zi, fixef(m_tmb)$zi, tolerance = tol_fixef)
})

test_that("ordbeta: prediction and simulation", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    family = ordbeta(),
    data = ordbeta_dat,
    start = ordbeta_start,
    map = ordbeta_map,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    family = ordbeta(),
    data = ordbeta_dat,
    start = ordbeta_start,
    map = ordbeta_map,
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

  set.seed(2003)
  sim_y <- simulate(m_rtmb, nsim = 1)[[1]]
  expect_equal(length(sim_y), nrow(ordbeta_dat))
  expect_true(all(sim_y >= 0 & sim_y <= 1))
})
