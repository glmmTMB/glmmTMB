## Test cases for the RTMB zero-truncated compois family
## Fit each model with RTMB and legacy TMB, then compare likelihoods,
## fixed effects, and covariance estimates where applicable.

context("RTMB truncated compois backend")

skip_if_not_installed("RTMB")

tol_logLik <- 1e-5
tol_fixef <- 1e-5
tol_varcorr <- 1e-4

cmpdat <- data.frame(
  f = factor(rep(c("a", "b"), 10)),
  y = c(15, 5, 20, 7, 19, 7, 19, 7, 19, 6,
        19, 10, 20, 8, 21, 8, 22, 7, 20, 8)
)

set.seed(1001)
truncated_compois_dat <- expand.grid(
  obs = seq_len(4),
  g = factor(seq_len(16))
)
truncated_compois_dat$x <- rnorm(nrow(truncated_compois_dat))
cond_effect <- rnorm(nlevels(truncated_compois_dat$g), sd = 0.25)
truncated_compois_mu <- exp(
  0.8 + 0.2 * truncated_compois_dat$x +
    cond_effect[truncated_compois_dat$g]
)
truncated_compois_dat$count <- glmmTMB:::rtruncated_compois2_rtmb(
  nrow(truncated_compois_dat),
  mean = truncated_compois_mu,
  nu = 1.2
)
truncated_compois_dat$off <- log(
  runif(nrow(truncated_compois_dat), 0.8, 1.2)
)
truncated_compois_dat$w <- runif(nrow(truncated_compois_dat), 0.5, 1.5)

truncated_compois_zi_dat <- truncated_compois_dat
set.seed(1002)
zero_rows <- sample(nrow(truncated_compois_zi_dat), 8)
truncated_compois_zi_dat$count[zero_rows] <- 0

test_that("truncated compois: RTMB density matches manual formula", {
  x <- 1:6
  mean <- seq(0.8, 3, length.out = length(x))
  nu <- seq(0.7, 1.5, length.out = length(x))
  log_nzprob <- glmmTMB:::log_nzprob_truncated_compois_rtmb(mean, nu)

  expect_equal(
    glmmTMB:::dtruncated_compois2_rtmb(x, mean, nu, log = TRUE),
    RTMB::dcompois2(x, mean = mean, nu = nu, log = TRUE) - log_nzprob,
    tolerance = 1e-12
  )
  expect_equal(
    glmmTMB:::dtruncated_compois2_rtmb(0, 2, 1.2, log = TRUE),
    -Inf
  )
})

test_that("truncated compois: fixed conditional effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ f,
    family = truncated_compois,
    data = cmpdat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ f,
    family = truncated_compois,
    data = cmpdat,
    se = FALSE
  )

  expect_equal(
    as.numeric(logLik(m_rtmb)),
    as.numeric(logLik(m_tmb)),
    tolerance = tol_logLik
  )
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = tol_fixef)
  expect_equal(sigma(m_rtmb), sigma(m_tmb), tolerance = tol_fixef)
})

test_that("truncated compois: offsets and weights", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x + offset(off),
    weights = w,
    family = truncated_compois,
    data = truncated_compois_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x + offset(off),
    weights = w,
    family = truncated_compois,
    data = truncated_compois_dat,
    se = FALSE
  )

  expect_equal(
    as.numeric(logLik(m_rtmb)),
    as.numeric(logLik(m_tmb)),
    tolerance = tol_logLik
  )
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = tol_fixef)
})

test_that("truncated compois: dispersion fixed effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x,
    dispformula = ~ x,
    family = truncated_compois,
    data = truncated_compois_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x,
    dispformula = ~ x,
    family = truncated_compois,
    data = truncated_compois_dat,
    se = FALSE
  )

  expect_equal(
    as.numeric(logLik(m_rtmb)),
    as.numeric(logLik(m_tmb)),
    tolerance = tol_logLik
  )
  expect_equal(fixef(m_rtmb)$disp, fixef(m_tmb)$disp, tolerance = tol_fixef)
})

test_that("truncated compois: conditional random intercept", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x + (1 | g),
    family = truncated_compois,
    data = truncated_compois_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x + (1 | g),
    family = truncated_compois,
    data = truncated_compois_dat,
    se = FALSE
  )

  expect_equal(
    as.numeric(logLik(m_rtmb)),
    as.numeric(logLik(m_tmb)),
    tolerance = tol_logLik
  )
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = tol_fixef)
  expect_equal(
    as.numeric(VarCorr(m_rtmb)$cond$g),
    as.numeric(VarCorr(m_tmb)$cond$g),
    tolerance = tol_varcorr
  )
})

test_that("truncated compois: zero-inflation fixed effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x,
    ziformula = ~ x,
    family = truncated_compois,
    data = truncated_compois_zi_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x,
    ziformula = ~ x,
    family = truncated_compois,
    data = truncated_compois_zi_dat,
    se = FALSE
  )

  expect_equal(
    as.numeric(logLik(m_rtmb)),
    as.numeric(logLik(m_tmb)),
    tolerance = tol_logLik
  )
  expect_equal(fixef(m_rtmb)$zi, fixef(m_tmb)$zi, tolerance = tol_fixef)
})

test_that("truncated compois: prediction values", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x,
    family = truncated_compois,
    data = truncated_compois_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ x,
    family = truncated_compois,
    data = truncated_compois_dat,
    se = FALSE
  )

  expect_equal(
    predict(m_rtmb, type = "response"),
    predict(m_tmb, type = "response"),
    tolerance = tol_fixef
  )
})

test_that("truncated compois: simulate works under RTMB backend", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ x + (1 | g),
    family = truncated_compois,
    data = truncated_compois_dat,
    se = FALSE
  )
  sim <- m_rtmb$obj$simulate(complete = TRUE)
  simulated <- sim$yobs

  expect_true(is.list(sim))
  expect_equal(length(simulated), nrow(truncated_compois_dat))
  expect_true(all(simulated >= 1))
  expect_true(all(abs(simulated - round(simulated)) < .Machine$double.eps))
})
