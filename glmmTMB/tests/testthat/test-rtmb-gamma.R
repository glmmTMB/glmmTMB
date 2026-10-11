## Test cases for the RTMB Gamma family
## Fit each model with RTMB and legacy TMB, then compare likelihoods,
## fixed effects, predictions, and covariance estimates where applicable.

context("RTMB Gamma backend")

skip_if_not_installed("RTMB")

tol_logLik <- 1e-5
tol_fixef <- 1e-5
tol_varcorr <- 1e-4

set.seed(1001)
gamma_dat <- expand.grid(
  obs = seq_len(5),
  g = factor(seq_len(30))
)
gamma_dat$x <- rnorm(nrow(gamma_dat))
gamma_dat$off <- log(runif(nrow(gamma_dat), 0.8, 1.2))
gamma_dat$w <- runif(nrow(gamma_dat), 0.5, 1.5)
gamma_effect <- rnorm(nlevels(gamma_dat$g), sd = 0.25)
gamma_mu <- exp(
  0.6 + 0.2 * gamma_dat$x + gamma_effect[gamma_dat$g]
)
gamma_shape <- 3
gamma_dat$y <- rgamma(
  nrow(gamma_dat),
  shape = gamma_shape,
  scale = gamma_mu / gamma_shape
)

set.seed(1002)
gamma_zi_dat <- gamma_dat
zi_effect <- rnorm(nlevels(gamma_zi_dat$g), sd = 0.3)
zi_prob <- plogis(-1.8 + zi_effect[gamma_zi_dat$g])
is_structural_zero <- rbinom(
  nrow(gamma_zi_dat),
  size = 1,
  prob = zi_prob
)
gamma_zi_dat$y <- ifelse(is_structural_zero == 1, 0, gamma_zi_dat$y)

test_that("Gamma: RTMB density matches C++ parameterization", {
  x <- seq(0.25, 2, length.out = 5)
  mu <- seq(0.5, 2.5, length.out = 5)
  phi <- seq(1.5, 4, length.out = 5)

  expect_equal(
    glmmTMB:::dgamma_rtmb(x, mean = mu, shape = phi, log = TRUE),
    dgamma(x, shape = phi, scale = mu / phi, log = TRUE),
    tolerance = 1e-12
  )
  expect_equal(
    glmmTMB:::dgamma_rtmb(0, mean = 1, shape = 2, log = TRUE),
    -Inf
  )
})

test_that("Gamma: fixed conditional effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    family = Gamma(link = "log"),
    data = gamma_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    family = Gamma(link = "log"),
    data = gamma_dat,
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

test_that("Gamma: offsets and weights", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x + offset(off),
    weights = w,
    family = Gamma(link = "log"),
    data = gamma_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x + offset(off),
    weights = w,
    family = Gamma(link = "log"),
    data = gamma_dat,
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

test_that("Gamma: dispersion fixed effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    dispformula = ~ x,
    family = Gamma(link = "log"),
    data = gamma_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    dispformula = ~ x,
    family = Gamma(link = "log"),
    data = gamma_dat,
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
  expect_equal(
    sigma(m_rtmb),
    sigma(m_tmb),
    tolerance = tol_fixef
  )
})

test_that("Gamma: conditional random intercept", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x + (1 | g),
    family = Gamma(link = "log"),
    data = gamma_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x + (1 | g),
    family = Gamma(link = "log"),
    data = gamma_dat,
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

test_that("Gamma: identity link", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ 1,
    family = Gamma(link = "identity"),
    data = gamma_dat,
    start = list(beta = mean(gamma_dat$y)),
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ 1,
    family = Gamma(link = "identity"),
    data = gamma_dat,
    start = list(beta = mean(gamma_dat$y)),
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

test_that("Gamma: zero-inflation fixed effects with ziGamma", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    ziformula = ~ x,
    family = ziGamma(link = "log"),
    data = gamma_zi_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    ziformula = ~ x,
    family = ziGamma(link = "log"),
    data = gamma_zi_dat,
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

test_that("Gamma: zero-inflation random effects with ziGamma", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    ziformula = ~ 0 + (1 | g),
    family = ziGamma(link = "log"),
    data = gamma_zi_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    ziformula = ~ 0 + (1 | g),
    family = ziGamma(link = "log"),
    data = gamma_zi_dat,
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

test_that("Gamma: prediction and simulation", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    family = Gamma(link = "log"),
    data = gamma_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    family = Gamma(link = "log"),
    data = gamma_dat,
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

  set.seed(1003)
  sim_y <- simulate(m_rtmb, nsim = 1)[[1]]
  expect_equal(length(sim_y), nrow(gamma_dat))
  expect_true(all(sim_y > 0))
})
