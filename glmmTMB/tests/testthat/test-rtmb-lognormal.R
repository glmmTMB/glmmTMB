## Test cases for the RTMB lognormal family
## Fit each model with RTMB and legacy TMB, then compare likelihoods,
## fixed effects, predictions, and covariance estimates where applicable.

context("RTMB lognormal backend")

skip_if_not_installed("RTMB")

tol_logLik <- 1e-5
tol_fixef <- 1e-5
tol_varcorr <- 1e-4

set.seed(1101)
lognormal_dat <- expand.grid(
  obs = seq_len(5),
  g = factor(seq_len(30))
)
lognormal_dat$x <- rnorm(nrow(lognormal_dat))
lognormal_dat$off <- log(runif(nrow(lognormal_dat), 0.8, 1.2))
lognormal_dat$w <- runif(nrow(lognormal_dat), 0.5, 1.5)
cond_effect <- rnorm(nlevels(lognormal_dat$g), sd = 0.25)
lognormal_mu <- exp(
  0.5 + 0.2 * lognormal_dat$x + cond_effect[lognormal_dat$g]
)
lognormal_sd <- 0.45
lognormal_logvar <- log1p((lognormal_sd / lognormal_mu)^2)
lognormal_dat$y <- rlnorm(
  nrow(lognormal_dat),
  meanlog = log(lognormal_mu) - lognormal_logvar / 2,
  sdlog = sqrt(lognormal_logvar)
)

set.seed(1102)
lognormal_zi_dat <- lognormal_dat
zi_effect <- rnorm(nlevels(lognormal_zi_dat$g), sd = 0.3)
zi_prob <- plogis(-1.8 + zi_effect[lognormal_zi_dat$g])
is_structural_zero <- rbinom(
  nrow(lognormal_zi_dat),
  size = 1,
  prob = zi_prob
)
lognormal_zi_dat$y <- ifelse(
  is_structural_zero == 1,
  0,
  lognormal_zi_dat$y
)

test_that("lognormal: RTMB density matches data-scale parameterization", {
  x <- seq(0.25, 2, length.out = 5)
  mu <- seq(0.5, 2.5, length.out = 5)
  sd <- seq(0.25, 0.75, length.out = 5)
  logvar <- log1p((sd / mu)^2)

  expect_equal(
    glmmTMB:::dlognormal_rtmb(x, mean = mu, sd = sd, log = TRUE),
    dlnorm(
      x,
      meanlog = log(mu) - logvar / 2,
      sdlog = sqrt(logvar),
      log = TRUE
    ),
    tolerance = 1e-12
  )
  expect_equal(
    glmmTMB:::dlognormal_rtmb(0, mean = 1, sd = 0.5, log = TRUE),
    -Inf
  )
})

test_that("lognormal: fixed conditional effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    family = lognormal(link = "log"),
    data = lognormal_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    family = lognormal(link = "log"),
    data = lognormal_dat,
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

test_that("lognormal: offsets and weights", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x + offset(off),
    weights = w,
    family = lognormal(link = "log"),
    data = lognormal_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x + offset(off),
    weights = w,
    family = lognormal(link = "log"),
    data = lognormal_dat,
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

test_that("lognormal: dispersion fixed effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    dispformula = ~ x,
    family = lognormal(link = "log"),
    data = lognormal_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    dispformula = ~ x,
    family = lognormal(link = "log"),
    data = lognormal_dat,
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

test_that("lognormal: conditional random intercept", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x + (1 | g),
    family = lognormal(link = "log"),
    data = lognormal_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x + (1 | g),
    family = lognormal(link = "log"),
    data = lognormal_dat,
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

test_that("lognormal: identity link", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ 1,
    family = lognormal(link = "identity"),
    data = lognormal_dat,
    start = list(beta = mean(lognormal_dat$y)),
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ 1,
    family = lognormal(link = "identity"),
    data = lognormal_dat,
    start = list(beta = mean(lognormal_dat$y)),
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

test_that("lognormal: zero-inflation fixed effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    ziformula = ~ x,
    family = lognormal(link = "log"),
    data = lognormal_zi_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    ziformula = ~ x,
    family = lognormal(link = "log"),
    data = lognormal_zi_dat,
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

test_that("lognormal: zero-inflation random effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    ziformula = ~ 0 + (1 | g),
    family = lognormal(link = "log"),
    data = lognormal_zi_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    ziformula = ~ 0 + (1 | g),
    family = lognormal(link = "log"),
    data = lognormal_zi_dat,
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

test_that("lognormal: prediction and simulation", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ x,
    family = lognormal(link = "log"),
    data = lognormal_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ x,
    family = lognormal(link = "log"),
    data = lognormal_dat,
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

  set.seed(1103)
  sim_y <- simulate(m_rtmb, nsim = 1)[[1]]
  expect_equal(length(sim_y), nrow(lognormal_dat))
  expect_true(all(sim_y > 0))
})
