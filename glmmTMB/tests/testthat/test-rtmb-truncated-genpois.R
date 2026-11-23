## Test cases for the RTMB zero-truncated generalized Poisson family
## Fit each model with the RTMB and legacy TMB backends, then compare
## likelihoods, fixed effects, and covariance estimates where applicable.

context("RTMB truncated genpois backend")

skip_if_not_installed("RTMB")

tol_logLik <- 1e-5
tol_fixef <- 1e-5
tol_varcorr <- 1e-4

set.seed(1501)
truncated_genpois_dat <- expand.grid(
  obs = seq_len(5),
  g = factor(seq_len(30))
)
truncated_genpois_dat$x <- rnorm(nrow(truncated_genpois_dat))
cond_effect <- rnorm(nlevels(truncated_genpois_dat$g), sd = 0.25)
mu <- exp(
  0.7 + 0.2 * truncated_genpois_dat$x +
    cond_effect[truncated_genpois_dat$g]
)
phi <- 1.4
theta <- mu / sqrt(phi)
lambda <- rep(1 - 1 / sqrt(phi), length.out = length(mu))
truncated_genpois_dat$count <- vapply(
  seq_along(mu),
  function(i) glmmTMB:::rtruncated_genpois_rtmb(theta[i], lambda[i]),
  numeric(1)
)
truncated_genpois_dat$off <- log(
  runif(nrow(truncated_genpois_dat), 0.8, 1.2)
)
truncated_genpois_dat$w <- runif(nrow(truncated_genpois_dat), 0.5, 1.5)

truncated_genpois_zi_dat <- truncated_genpois_dat
zi_prob <- plogis(-1.6 + 0.4 * truncated_genpois_zi_dat$x)
zi_zero <- rbinom(nrow(truncated_genpois_zi_dat), size = 1, prob = zi_prob)
truncated_genpois_zi_dat$count <- ifelse(
  zi_zero == 1,
  0,
  truncated_genpois_zi_dat$count
)

expect_truncated_genpois_equal <- function(m_rtmb, m_tmb,
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

test_that("truncated genpois: RTMB density matches translated formula", {
  x <- 1:8
  theta <- 1.4
  lambda <- 0.2
  log_nzprob <- log(1 - exp(-theta))
  target <- glmmTMB::dgenpois(
    x,
    lambda1 = theta,
    lambda2 = lambda,
    log = TRUE
  ) - log_nzprob

  expect_equal(
    glmmTMB:::dtruncated_genpois_rtmb(
      x,
      theta = theta,
      lambda = lambda,
      log = TRUE
    ),
    target,
    tolerance = 1e-12
  )
  expect_equal(
    glmmTMB:::dtruncated_genpois_rtmb(
      0,
      theta = theta,
      lambda = lambda,
      log = TRUE
    ),
    -Inf
  )
})

test_that("truncated genpois: fixed conditional effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x,
    family = truncated_genpois,
    data = truncated_genpois_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x,
    family = truncated_genpois,
    data = truncated_genpois_dat,
    se = FALSE
  )

  expect_truncated_genpois_equal(m_rtmb, m_tmb)
})

test_that("truncated genpois: offsets and weights", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x + offset(off),
    weights = w,
    family = truncated_genpois,
    data = truncated_genpois_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x + offset(off),
    weights = w,
    family = truncated_genpois,
    data = truncated_genpois_dat,
    se = FALSE
  )

  expect_truncated_genpois_equal(m_rtmb, m_tmb)
})

test_that("truncated genpois: dispersion fixed effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x,
    dispformula = ~ x,
    family = truncated_genpois,
    data = truncated_genpois_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x,
    dispformula = ~ x,
    family = truncated_genpois,
    data = truncated_genpois_dat,
    se = FALSE
  )

  expect_truncated_genpois_equal(m_rtmb, m_tmb)
})

test_that("truncated genpois: conditional random intercept", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x + (1 | g),
    family = truncated_genpois,
    data = truncated_genpois_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x + (1 | g),
    family = truncated_genpois,
    data = truncated_genpois_dat,
    se = FALSE
  )

  expect_truncated_genpois_equal(m_rtmb, m_tmb, check_varcorr = TRUE)
})

test_that("truncated genpois: fixed zero inflation", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x,
    ziformula = ~ x,
    family = truncated_genpois,
    data = truncated_genpois_zi_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x,
    ziformula = ~ x,
    family = truncated_genpois,
    data = truncated_genpois_zi_dat,
    se = FALSE
  )

  expect_truncated_genpois_equal(m_rtmb, m_tmb)
})

test_that("truncated genpois: prediction with standard errors", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x,
    family = truncated_genpois,
    data = truncated_genpois_dat
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x,
    family = truncated_genpois,
    data = truncated_genpois_dat
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

test_that("truncated genpois: simulation is strictly positive", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x,
    family = truncated_genpois,
    data = truncated_genpois_dat,
    se = FALSE
  )

  set.seed(1502)
  sim <- simulate(m_rtmb, nsim = 1)
  simulated <- sim[[1]]

  expect_length(simulated, nrow(truncated_genpois_dat))
  expect_true(all(simulated >= 1))
  expect_true(all(simulated == floor(simulated)))
})
