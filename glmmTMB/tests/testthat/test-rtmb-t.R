## Test cases for the RTMB t family
## Fit each model with RTMB and legacy TMB, then compare likelihoods,
## fixed effects, predictions, and covariance estimates where applicable.

context("RTMB t backend")

skip_if_not_installed("RTMB")

tol_logLik <- 1e-5
tol_fixef <- 1e-5
tol_varcorr <- 1e-4

set.seed(1201)
t_dat <- expand.grid(
  obs = seq_len(5),
  g = factor(seq_len(30))
)
t_dat$x <- rnorm(nrow(t_dat))
t_dat$off <- rnorm(nrow(t_dat), sd = 0.1)
t_dat$w <- runif(nrow(t_dat), 0.5, 1.5)
t_effect <- rnorm(nlevels(t_dat$g), sd = 0.25)
t_mu <- 0.4 + 0.3 * t_dat$x + t_effect[t_dat$g]
t_scale <- 0.7
t_df <- 6
t_dat$y <- t_mu + t_scale * rt(nrow(t_dat), df = t_df)

set.seed(1202)
t_zi_dat <- t_dat
zi_effect <- rnorm(nlevels(t_zi_dat$g), sd = 0.3)
zi_prob <- plogis(-1.8 + zi_effect[t_zi_dat$g])
is_structural_zero <- rbinom(nrow(t_zi_dat), size = 1, prob = zi_prob)
t_zi_dat$y <- ifelse(is_structural_zero == 1, 0, t_zi_dat$y)

t_start <- list(psi = log(t_df))
t_map <- list(psi = factor(NA))

expect_t_matches <- function(m_rtmb, m_tmb, check_varcorr = FALSE) {
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

test_that("t: RTMB density matches C++ parameterization", {
  x <- seq(-1, 1, length.out = 5)
  mu <- seq(-0.2, 0.2, length.out = 5)
  scale <- seq(0.5, 1.5, length.out = 5)
  df <- 6

  expect_equal(
    glmmTMB:::dt_rtmb(x, mean = mu, scale = scale, df = df, log = TRUE),
    dt((x - mu) / scale, df = df, log = TRUE) - log(scale),
    tolerance = 1e-12
  )
})

test_that("t: fixed conditional effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    family = t_family,
    data = t_dat,
    start = t_start,
    map = t_map,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    family = t_family,
    data = t_dat,
    start = t_start,
    map = t_map,
    se = FALSE
  )

  expect_t_matches(m_rtmb, m_tmb)
})

test_that("t: offsets and weights", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x + offset(off),
    weights = w,
    family = t_family,
    data = t_dat,
    start = t_start,
    map = t_map,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x + offset(off),
    weights = w,
    family = t_family,
    data = t_dat,
    start = t_start,
    map = t_map,
    se = FALSE
  )

  expect_t_matches(m_rtmb, m_tmb)
})

test_that("t: dispersion fixed effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    dispformula = ~ x,
    family = t_family,
    data = t_dat,
    start = t_start,
    map = t_map,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    dispformula = ~ x,
    family = t_family,
    data = t_dat,
    start = t_start,
    map = t_map,
    se = FALSE
  )

  expect_t_matches(m_rtmb, m_tmb)
  expect_equal(
    fixef(m_rtmb)$disp,
    fixef(m_tmb)$disp,
    tolerance = tol_fixef
  )
  expect_equal(
    sigma(m_rtmb),
    sigma(m_tmb),
    tolerance = tol_fixef
  )
})

test_that("t: conditional random intercept", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x + (1 | g),
    family = t_family,
    data = t_dat,
    start = t_start,
    map = t_map,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x + (1 | g),
    family = t_family,
    data = t_dat,
    start = t_start,
    map = t_map,
    se = FALSE
  )

  expect_t_matches(m_rtmb, m_tmb, check_varcorr = TRUE)
})

test_that("t: estimated degrees of freedom parameter", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    family = t_family,
    data = t_dat,
    start = t_start,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    family = t_family,
    data = t_dat,
    start = t_start,
    se = FALSE
  )

  expect_t_matches(m_rtmb, m_tmb)
  expect_equal(
    m_rtmb$obj$env$last.par.best["psi"],
    m_tmb$obj$env$last.par.best["psi"],
    tolerance = tol_fixef
  )
})

test_that("t: zero-inflation fixed effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    ziformula = ~ x,
    family = t_family,
    data = t_zi_dat,
    start = t_start,
    map = t_map,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    ziformula = ~ x,
    family = t_family,
    data = t_zi_dat,
    start = t_start,
    map = t_map,
    se = FALSE
  )

  expect_t_matches(m_rtmb, m_tmb)
  expect_equal(
    fixef(m_rtmb)$zi,
    fixef(m_tmb)$zi,
    tolerance = tol_fixef
  )
})

test_that("t: prediction and simulation", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    family = t_family,
    data = t_dat,
    start = t_start,
    map = t_map,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    family = t_family,
    data = t_dat,
    start = t_start,
    map = t_map,
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

  set.seed(1203)
  sim_y <- simulate(m_rtmb, nsim = 1)[[1]]
  expect_equal(length(sim_y), nrow(t_dat))
  expect_true(all(is.finite(sim_y)))
})
