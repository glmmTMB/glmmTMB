## Test cases for the RTMB zero-truncated Poisson family
## Fit each model with RTMB and legacy TMB, then compare likelihoods,
## fixed effects, and covariance estimates where applicable.

context("RTMB truncated Poisson backend")

skip_if_not_installed("RTMB")

data("Salamanders", package = "glmmTMB")

tol_logLik <- 1e-5
tol_fixef <- 1e-5
tol_varcorr <- 1e-4

positive_salamanders <- subset(Salamanders, count > 0)

test_that("truncated Poisson: RTMB density matches package density", {
  x <- 1:8
  lambda <- seq(0.25, 4, length.out = length(x))

  expect_equal(
    glmmTMB:::dtruncated_poisson_rtmb(x, lambda, log = TRUE),
    glmmTMB::dtruncated_poisson(x, lambda, log = TRUE),
    tolerance = 1e-12
  )
  expect_equal(
    glmmTMB:::dtruncated_poisson_rtmb(x, lambda),
    glmmTMB::dtruncated_poisson(x, lambda),
    tolerance = 1e-12
  )
  expect_equal(
    glmmTMB:::dtruncated_poisson_rtmb(0, 1, log = TRUE),
    -Inf
  )
})

test_that("truncated Poisson: fixed conditional effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ mined,
    family = truncated_poisson,
    data = positive_salamanders,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ mined,
    family = truncated_poisson,
    data = positive_salamanders,
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

test_that("truncated Poisson: conditional offset", {
  set.seed(3001)
  offset_data <- transform(
    positive_salamanders,
    log_exposure = log(runif(nrow(positive_salamanders), 0.5, 2))
  )

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ mined + offset(log_exposure),
    family = truncated_poisson,
    data = offset_data,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ mined + offset(log_exposure),
    family = truncated_poisson,
    data = offset_data,
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

test_that("truncated Poisson: observation weights", {
  weighted_data <- transform(
    positive_salamanders,
    w = rep(c(0.5, 1, 2), length.out = nrow(positive_salamanders))
  )

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ mined,
    family = truncated_poisson,
    weights = w,
    data = weighted_data,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ mined,
    family = truncated_poisson,
    weights = w,
    data = weighted_data,
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

test_that("truncated Poisson: identity link", {
  identity_start <- mean(positive_salamanders$count)

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ 1,
    family = truncated_poisson(link = "identity"),
    data = positive_salamanders,
    start = list(beta = identity_start),
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ 1,
    family = truncated_poisson(link = "identity"),
    data = positive_salamanders,
    start = list(beta = identity_start),
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

test_that("truncated Poisson: square-root link", {
  sqrt_start <- sqrt(mean(positive_salamanders$count))

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ 1,
    family = truncated_poisson(link = "sqrt"),
    data = positive_salamanders,
    start = list(beta = sqrt_start),
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ 1,
    family = truncated_poisson(link = "sqrt"),
    data = positive_salamanders,
    start = list(beta = sqrt_start),
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

test_that("truncated Poisson: conditional random intercept", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ mined + (1 | site),
    family = truncated_poisson,
    data = positive_salamanders,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ mined + (1 | site),
    family = truncated_poisson,
    data = positive_salamanders,
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
    VarCorr(m_rtmb),
    VarCorr(m_tmb),
    tolerance = tol_varcorr
  )
})

test_that("truncated Poisson: sparse conditional matrix", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ mined + spp,
    family = truncated_poisson,
    data = positive_salamanders,
    sparseX = c(cond = TRUE),
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ mined + spp,
    family = truncated_poisson,
    data = positive_salamanders,
    sparseX = c(cond = TRUE),
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

test_that("truncated Poisson: fixed zero inflation", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ mined,
    ziformula = ~ mined,
    family = truncated_poisson,
    data = Salamanders,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ mined,
    ziformula = ~ mined,
    family = truncated_poisson,
    data = Salamanders,
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
    fixef(m_rtmb)$zi,
    fixef(m_tmb)$zi,
    tolerance = tol_fixef
  )
})

test_that("truncated Poisson: zero-inflation random intercept", {
  thetazi <- log(0.5)
  theta_map <- factor(NA)

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ mined,
    ziformula = ~ mined + (1 | site),
    family = truncated_poisson,
    data = Salamanders,
    start = list(thetazi = thetazi),
    map = list(thetazi = theta_map),
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ mined,
    ziformula = ~ mined + (1 | site),
    family = truncated_poisson,
    data = Salamanders,
    start = list(thetazi = thetazi),
    map = list(thetazi = theta_map),
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
  expect_equal(
    VarCorr(m_rtmb),
    VarCorr(m_tmb),
    tolerance = tol_varcorr
  )
})

test_that("truncated Poisson: random-only zero inflation", {
  thetazi <- log(0.5)
  theta_map <- factor(NA)

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ mined,
    ziformula = ~ 0 + (1 | site),
    family = truncated_poisson,
    data = Salamanders,
    start = list(thetazi = thetazi),
    map = list(thetazi = theta_map),
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ mined,
    ziformula = ~ 0 + (1 | site),
    family = truncated_poisson,
    data = Salamanders,
    start = list(thetazi = thetazi),
    map = list(thetazi = theta_map),
    se = FALSE
  )

  expect_equal(
    as.numeric(logLik(m_rtmb)),
    as.numeric(logLik(m_tmb)),
    tolerance = tol_logLik
  )
  expect_equal(
    VarCorr(m_rtmb),
    VarCorr(m_tmb),
    tolerance = tol_varcorr
  )
})

test_that("truncated Poisson: simulation is strictly positive", {
  local_useRTMB(TRUE)
  model <- glmmTMB(
    count ~ mined + (1 | site),
    family = truncated_poisson,
    data = positive_salamanders,
    se = FALSE
  )

  set.seed(3002)
  simulated <- model$obj$simulate(complete = TRUE)$yobs

  expect_length(simulated, nrow(positive_salamanders))
  expect_true(all(is.finite(simulated)))
  expect_true(all(simulated > 0))
  expect_true(all(simulated == floor(simulated)))
})

test_that("truncated Poisson: simulation remains positive for tiny means", {
  local_useRTMB(TRUE)
  model <- glmmTMB(
    count ~ 1,
    family = truncated_poisson,
    data = positive_salamanders,
    start = list(beta = log(1e-8)),
    map = list(beta = factor(NA)),
    se = FALSE
  )

  set.seed(3003)
  simulated <- model$obj$simulate(complete = TRUE)$yobs

  expect_true(all(is.finite(simulated)))
  expect_true(all(simulated > 0))
})

test_that("truncated Poisson: ZI simulation contains zeros and positives", {
  fixed_map <- factor(NA)

  local_useRTMB(TRUE)
  model <- glmmTMB(
    count ~ 1,
    ziformula = ~ 1,
    family = truncated_poisson,
    data = Salamanders,
    start = list(beta = log(10), betazi = qlogis(0.5)),
    map = list(beta = fixed_map, betazi = fixed_map),
    se = FALSE
  )

  set.seed(3004)
  simulated <- model$obj$simulate(complete = TRUE)$yobs

  expect_length(simulated, nrow(Salamanders))
  expect_true(all(is.finite(simulated)))
  expect_true(any(simulated == 0))
  expect_true(any(simulated > 0))
  expect_true(all(simulated == floor(simulated)))
})
