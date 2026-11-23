## Test cases for the RTMB Bell family
## Fit each model with RTMB and legacy TMB, then compare likelihoods,
## parameters, predictions, and covariance estimates where applicable.

context("RTMB Bell backend")

skip_if_not_installed("RTMB")
skip_if_not_installed("gsl")

tol_logLik <- 1e-5
tol_fixef <- 1e-5
tol_varcorr <- 1e-4

rbell_test <- function(theta) {
  vapply(theta, function(theta_i) {
    n_compound <- rpois(1L, expm1(theta_i))
    if (n_compound == 0L) {
      return(0)
    }
    sum(vapply(seq_len(n_compound), function(i) {
      repeat {
        value <- rpois(1L, theta_i)
        if (value > 0L) {
          return(value)
        }
      }
    }, numeric(1)))
  }, numeric(1))
}

set.seed(2201)
bell_dat <- expand.grid(
  obs = seq_len(8),
  g = factor(seq_len(40))
)
bell_dat$x <- rnorm(nrow(bell_dat))
bell_dat$off <- runif(nrow(bell_dat), -0.15, 0.15)
bell_dat$w <- runif(nrow(bell_dat), 0.5, 1.5)
cond_effect <- rnorm(nlevels(bell_dat$g), sd = 0.2)
bell_mu <- exp(0.7 + 0.25 * bell_dat$x + cond_effect[bell_dat$g])
bell_theta <- gsl::lambert_W0(bell_mu)
bell_dat$y <- rbell_test(bell_theta)

set.seed(2202)
bell_zi_dat <- bell_dat
zi_effect <- rnorm(nlevels(bell_zi_dat$g), sd = 0.7)
zi_prob <- plogis(-1.3 + zi_effect[bell_zi_dat$g])
is_structural_zero <- rbinom(nrow(bell_zi_dat), 1L, zi_prob)
bell_zi_dat$y[is_structural_zero == 1L] <- 0

expect_bell_matches <- function(m_rtmb, m_tmb, check_varcorr = FALSE) {
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
  if (check_varcorr) {
    expect_equal(
      as.numeric(VarCorr(m_rtmb)$cond$g),
      as.numeric(VarCorr(m_tmb)$cond$g),
      tolerance = tol_varcorr
    )
  }
}

test_that("bell: RTMB density matches package density", {
  x <- 0:8
  theta <- seq(0.2, 1.2, length.out = length(x))
  mean <- theta * exp(theta)

  expect_equal(
    glmmTMB:::dbell_rtmb(x, mean = mean, log = TRUE),
    glmmTMB::dbell(x, theta = theta, log = TRUE),
    tolerance = 1e-12
  )
})

test_that("bell: fixed conditional effects, offsets, and weights", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x + offset(off),
    weights = w,
    family = bell(),
    data = bell_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x + offset(off),
    weights = w,
    family = bell(),
    data = bell_dat,
    se = FALSE
  )

  expect_bell_matches(m_rtmb, m_tmb)
})

test_that("bell: conditional random intercept", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x + (1 | g),
    family = bell(),
    data = bell_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x + (1 | g),
    family = bell(),
    data = bell_dat,
    se = FALSE
  )

  expect_bell_matches(m_rtmb, m_tmb, check_varcorr = TRUE)
})

test_that("bell: zero-inflation fixed effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    ziformula = ~ x,
    family = bell(),
    data = bell_zi_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    ziformula = ~ x,
    family = bell(),
    data = bell_zi_dat,
    se = FALSE
  )

  expect_bell_matches(m_rtmb, m_tmb)
})

test_that("bell: zero-inflation random effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    ziformula = ~ 1 + (1 | g),
    family = bell(),
    data = bell_zi_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    ziformula = ~ 1 + (1 | g),
    family = bell(),
    data = bell_zi_dat,
    se = FALSE
  )

  expect_bell_matches(m_rtmb, m_tmb)
  expect_equal(
    as.numeric(VarCorr(m_rtmb)$zi$g),
    as.numeric(VarCorr(m_tmb)$zi$g),
    tolerance = tol_varcorr
  )
})

test_that("bell: prediction and simulation", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    family = bell(),
    data = bell_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    family = bell(),
    data = bell_dat,
    se = FALSE
  )

  expect_equal(
    predict(m_rtmb, type = "response"),
    predict(m_tmb, type = "response"),
    tolerance = tol_fixef
  )

  set.seed(2203)
  simulated <- simulate(m_rtmb, nsim = 1)[[1]]
  expect_equal(length(simulated), nrow(bell_dat))
  expect_true(all(simulated >= 0))
  expect_true(all(simulated == floor(simulated)))
})
