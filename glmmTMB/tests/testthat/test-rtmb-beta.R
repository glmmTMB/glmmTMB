## Test cases for the RTMB beta family
## Fit each model with the RTMB and legacy TMB backends, then compare
## likelihoods, fixed effects, predictions, and covariance estimates where
## applicable.

context("RTMB beta backend")

skip_if_not_installed("RTMB")

tol_logLik <- 1e-5
tol_fixef <- 1e-5
tol_varcorr <- 1e-4

set.seed(1801)
beta_dat <- expand.grid(
  obs = seq_len(5),
  g = factor(seq_len(30))
)
beta_dat$x <- rnorm(nrow(beta_dat))
beta_dat$off <- runif(nrow(beta_dat), -0.2, 0.2)
beta_dat$w <- runif(nrow(beta_dat), 0.5, 1.5)
beta_effect <- rnorm(nlevels(beta_dat$g), sd = 0.25)
beta_mu <- plogis(
  -0.2 + 0.35 * beta_dat$x + beta_effect[beta_dat$g]
)
beta_phi <- 8
beta_dat$y <- rbeta(
  nrow(beta_dat),
  shape1 = beta_mu * beta_phi,
  shape2 = (1 - beta_mu) * beta_phi
)

set.seed(1802)
beta_zi_dat <- beta_dat
zi_effect <- rnorm(nlevels(beta_zi_dat$g), sd = 0.35)
zi_prob <- plogis(-1.8 + zi_effect[beta_zi_dat$g])
is_structural_zero <- rbinom(
  nrow(beta_zi_dat),
  size = 1,
  prob = zi_prob
)
beta_zi_dat$y <- ifelse(is_structural_zero == 1, 0, beta_zi_dat$y)

test_that("beta: RTMB density matches C++ parameterization", {
  x <- seq(0.1, 0.9, length.out = 5)
  mu <- seq(0.2, 0.8, length.out = 5)
  phi <- seq(3, 7, length.out = 5)

  expect_equal(
    glmmTMB:::dbeta_rtmb(x, mean = mu, phi = phi, log = TRUE),
    dbeta(x, shape1 = mu * phi, shape2 = (1 - mu) * phi, log = TRUE),
    tolerance = 1e-12
  )
  expect_equal(
    glmmTMB:::dbeta_rtmb(0, mean = 0.4, phi = 5, log = TRUE),
    -Inf
  )
})

test_that("beta: fixed conditional effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    family = beta_family(link = "logit"),
    data = beta_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    family = beta_family(link = "logit"),
    data = beta_dat,
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

test_that("beta: offsets and weights", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x + offset(off),
    weights = w,
    family = beta_family(link = "logit"),
    data = beta_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x + offset(off),
    weights = w,
    family = beta_family(link = "logit"),
    data = beta_dat,
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

test_that("beta: dispersion fixed effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    dispformula = ~ x,
    family = beta_family(link = "logit"),
    data = beta_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    dispformula = ~ x,
    family = beta_family(link = "logit"),
    data = beta_dat,
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

test_that("beta: conditional random intercept", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x + (1 | g),
    family = beta_family(link = "logit"),
    data = beta_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x + (1 | g),
    family = beta_family(link = "logit"),
    data = beta_dat,
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

test_that("beta: probit link", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    family = beta_family(link = "probit"),
    data = beta_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    family = beta_family(link = "probit"),
    data = beta_dat,
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

test_that("beta: zero-inflation fixed effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    ziformula = ~ x,
    family = beta_family(link = "logit"),
    data = beta_zi_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    ziformula = ~ x,
    family = beta_family(link = "logit"),
    data = beta_zi_dat,
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

test_that("beta: zero-inflation random effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    ziformula = ~ 0 + (1 | g),
    family = beta_family(link = "logit"),
    data = beta_zi_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    ziformula = ~ 0 + (1 | g),
    family = beta_family(link = "logit"),
    data = beta_zi_dat,
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

test_that("beta: prediction and simulation", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    family = beta_family(link = "logit"),
    data = beta_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    family = beta_family(link = "logit"),
    data = beta_dat,
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

  set.seed(1803)
  sim_y <- simulate(m_rtmb, nsim = 1)[[1]]
  expect_equal(length(sim_y), nrow(beta_dat))
  expect_true(all(sim_y > 0 & sim_y < 1))
})
