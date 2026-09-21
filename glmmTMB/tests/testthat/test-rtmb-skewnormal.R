## Test cases for the RTMB skewnormal family
## Fit each model with RTMB and legacy TMB, then compare likelihoods,
## fixed effects, predictions, and covariance estimates where applicable.

context("RTMB skewnormal backend")

skip_if_not_installed("RTMB")

tol_logLik <- 1e-5
tol_fixef <- 1e-5
tol_varcorr <- 1e-4

set.seed(1301)
skewnormal_dat <- expand.grid(
  obs = seq_len(5),
  g = factor(seq_len(30))
)
skewnormal_dat$x <- rnorm(nrow(skewnormal_dat))
skewnormal_dat$off <- rnorm(nrow(skewnormal_dat), sd = 0.1)
skewnormal_dat$w <- runif(nrow(skewnormal_dat), 0.5, 1.5)
cond_effect <- rnorm(nlevels(skewnormal_dat$g), sd = 0.25)
skew_mu <- 0.4 + 0.3 * skewnormal_dat$x + cond_effect[skewnormal_dat$g]
skew_sd <- 0.8
skew_alpha <- -2
delta <- skew_alpha / sqrt(1 + skew_alpha^2)
omega <- skew_sd / sqrt(1 - 2 / pi * delta^2)
xi <- skew_mu - omega * delta * sqrt(2 / pi)
skewnormal_dat$y <- xi + omega * (
  delta * abs(rnorm(nrow(skewnormal_dat))) +
    sqrt(1 - delta^2) * rnorm(nrow(skewnormal_dat))
)

set.seed(1302)
skewnormal_zi_dat <- skewnormal_dat
zi_effect <- rnorm(nlevels(skewnormal_zi_dat$g), sd = 0.3)
zi_prob <- plogis(-1.8 + zi_effect[skewnormal_zi_dat$g])
is_structural_zero <- rbinom(
  nrow(skewnormal_zi_dat),
  size = 1,
  prob = zi_prob
)
skewnormal_zi_dat$y <- ifelse(
  is_structural_zero == 1,
  0,
  skewnormal_zi_dat$y
)

skew_start <- list(psi = skew_alpha)
skew_map <- list(psi = factor(NA))

expect_skewnormal_matches <- function(m_rtmb, m_tmb,
                                      check_varcorr = FALSE) {
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

test_that("skewnormal: RTMB density matches translated formula", {
  x <- seq(-1, 1, length.out = 5)
  mu <- seq(-0.2, 0.2, length.out = 5)
  sd <- seq(0.5, 1.5, length.out = 5)
  alpha <- -2
  delta <- alpha / sqrt(1 + alpha^2)
  omega <- sd / sqrt(1 - 2 / pi * delta^2)
  xi <- mu - omega * delta * sqrt(2 / pi)
  z <- (x - xi) / omega

  expect_equal(
    glmmTMB:::dskewnormal_rtmb(
      x,
      mean = mu,
      sd = sd,
      alpha = alpha,
      log = TRUE
    ),
    log(2) - log(omega) + dnorm(z, log = TRUE) + log(pnorm(alpha * z)),
    tolerance = 1e-12
  )
})

test_that("skewnormal: fixed conditional effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    family = skewnormal,
    data = skewnormal_dat,
    start = skew_start,
    map = skew_map,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    family = skewnormal,
    data = skewnormal_dat,
    start = skew_start,
    map = skew_map,
    se = FALSE
  )

  expect_skewnormal_matches(m_rtmb, m_tmb)
})

test_that("skewnormal: offsets and weights", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x + offset(off),
    weights = w,
    family = skewnormal,
    data = skewnormal_dat,
    start = skew_start,
    map = skew_map,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x + offset(off),
    weights = w,
    family = skewnormal,
    data = skewnormal_dat,
    start = skew_start,
    map = skew_map,
    se = FALSE
  )

  expect_skewnormal_matches(m_rtmb, m_tmb)
})

test_that("skewnormal: dispersion fixed effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    dispformula = ~ x,
    family = skewnormal,
    data = skewnormal_dat,
    start = skew_start,
    map = skew_map,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    dispformula = ~ x,
    family = skewnormal,
    data = skewnormal_dat,
    start = skew_start,
    map = skew_map,
    se = FALSE
  )

  expect_skewnormal_matches(m_rtmb, m_tmb)
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

test_that("skewnormal: conditional random intercept", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x + (1 | g),
    family = skewnormal,
    data = skewnormal_dat,
    start = skew_start,
    map = skew_map,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x + (1 | g),
    family = skewnormal,
    data = skewnormal_dat,
    start = skew_start,
    map = skew_map,
    se = FALSE
  )

  expect_skewnormal_matches(m_rtmb, m_tmb, check_varcorr = TRUE)
})

test_that("skewnormal: estimated shape parameter", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    family = skewnormal,
    data = skewnormal_dat,
    start = skew_start,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    family = skewnormal,
    data = skewnormal_dat,
    start = skew_start,
    se = FALSE
  )

  expect_skewnormal_matches(m_rtmb, m_tmb)
  expect_equal(
    m_rtmb$obj$env$last.par.best["psi"],
    m_tmb$obj$env$last.par.best["psi"],
    tolerance = tol_fixef
  )
})

test_that("skewnormal: zero-inflation fixed effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    ziformula = ~ x,
    family = skewnormal,
    data = skewnormal_zi_dat,
    start = skew_start,
    map = skew_map,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    ziformula = ~ x,
    family = skewnormal,
    data = skewnormal_zi_dat,
    start = skew_start,
    map = skew_map,
    se = FALSE
  )

  expect_skewnormal_matches(m_rtmb, m_tmb)
  expect_equal(
    fixef(m_rtmb)$zi,
    fixef(m_tmb)$zi,
    tolerance = tol_fixef
  )
})

test_that("skewnormal: prediction and simulation", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    family = skewnormal,
    data = skewnormal_dat,
    start = skew_start,
    map = skew_map,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    family = skewnormal,
    data = skewnormal_dat,
    start = skew_start,
    map = skew_map,
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

  set.seed(1303)
  sim_y <- simulate(m_rtmb, nsim = 1)[[1]]
  expect_equal(length(sim_y), nrow(skewnormal_dat))
  expect_true(all(is.finite(sim_y)))
})
