## Test cases for the RTMB Poisson family
## Fit each model with the RTMB and legacy TMB backends, then compare
## likelihoods, fixed effects, and covariance estimates where applicable.

data("Salamanders", package = "glmmTMB")

tol_logLik <- 1e-5
tol_fixef <- 1e-5
tol_varcorr <- 1e-4

set.seed(2001)
poisson_panel <- expand.grid(
  time = 0:3,
  group = factor(seq_len(30))
)
poisson_panel$time_fac <- factor(poisson_panel$time)
poisson_panel$time_num <- glmmTMB::numFactor(poisson_panel$time)
panel_effect <- rnorm(nlevels(poisson_panel$group), sd = 0.4)
poisson_panel$count <- rpois(
  nrow(poisson_panel),
  exp(
    0.5 + 0.15 * poisson_panel$time +
      panel_effect[poisson_panel$group]
  )
)

test_that("poisson: fixed conditional effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ mined,
    family = poisson,
    data = Salamanders,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ mined,
    family = poisson,
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
})

test_that("poisson: multi-level factor fixed effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ spp,
    family = poisson,
    data = Salamanders,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ spp,
    family = poisson,
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
})

test_that("poisson: conditional offset", {
  offset_data <- transform(
    Salamanders,
    log_exposure = log(runif(nrow(Salamanders), 0.5, 2))
  )

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ mined + offset(log_exposure),
    family = poisson,
    data = offset_data,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ mined + offset(log_exposure),
    family = poisson,
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

test_that("poisson: observation weights", {
  weighted_data <- transform(
    Salamanders,
    w = rep(c(0.5, 1, 2), length.out = nrow(Salamanders))
  )

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ mined,
    family = poisson,
    weights = w,
    data = weighted_data,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ mined,
    family = poisson,
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

test_that("poisson: missing responses", {
  missing_data <- Salamanders
  missing_data$count[c(3, 21, 55)] <- NA

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ mined,
    family = poisson,
    data = missing_data,
    na.action = na.exclude,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ mined,
    family = poisson,
    data = missing_data,
    na.action = na.exclude,
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

test_that("poisson: sparse conditional fixed-effects matrix", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ mined + spp,
    family = poisson,
    data = Salamanders,
    sparseX = c(cond = TRUE),
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ mined + spp,
    family = poisson,
    data = Salamanders,
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

test_that("poisson: explicit log link", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ mined,
    family = poisson(link = "log"),
    data = Salamanders,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ mined,
    family = poisson(link = "log"),
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
})

test_that("poisson: identity link", {
  identity_start <- mean(Salamanders$count)

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ 1,
    family = poisson(link = "identity"),
    data = Salamanders,
    start = list(beta = identity_start),
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ 1,
    family = poisson(link = "identity"),
    data = Salamanders,
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

test_that("poisson: square-root link", {
  sqrt_start <- sqrt(mean(Salamanders$count))

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ 1,
    family = poisson(link = "sqrt"),
    data = Salamanders,
    start = list(beta = sqrt_start),
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ 1,
    family = poisson(link = "sqrt"),
    data = Salamanders,
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

test_that("poisson: conditional random intercept", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ mined + (1 | site),
    family = poisson,
    data = Salamanders,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ mined + (1 | site),
    family = poisson,
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
    VarCorr(m_rtmb),
    VarCorr(m_tmb),
    tolerance = tol_varcorr
  )
})

test_that("poisson: estfun retapes with updated RTMB data", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ DOP + (1 | sample),
    family = poisson,
    data = Salamanders
  )
  result_rtmb <- expect_silent(sandwich::estfun(m_rtmb))

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ DOP + (1 | sample),
    family = poisson,
    data = Salamanders
  )
  result_tmb <- expect_silent(sandwich::estfun(m_tmb))

  expect_gt(max(abs(result_rtmb)), 1)
  expect_equal(result_rtmb, result_tmb, tolerance = 1e-3)
})

test_that("poisson: correlated conditional random slope", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ time + (time | group),
    family = poisson,
    data = poisson_panel,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ time + (time | group),
    family = poisson,
    data = poisson_panel,
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

test_that("poisson: multiple conditional random-effect terms", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ mined + (1 | site) + (1 | spp),
    family = poisson,
    data = Salamanders,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ mined + (1 | site) + (1 | spp),
    family = poisson,
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
    VarCorr(m_rtmb),
    VarCorr(m_tmb),
    tolerance = tol_varcorr
  )
})

test_that("poisson: diagonal conditional covariance", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ time + diag(time | group),
    family = poisson,
    data = poisson_panel,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ time + diag(time | group),
    family = poisson,
    data = poisson_panel,
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

test_that("poisson: homogeneous diagonal conditional covariance", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ time + homdiag(time | group),
    family = poisson,
    data = poisson_panel,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ time + homdiag(time | group),
    family = poisson,
    data = poisson_panel,
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

test_that("poisson: fixed zero-inflation intercept", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ mined,
    ziformula = ~ 1,
    family = poisson,
    data = Salamanders,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ mined,
    ziformula = ~ 1,
    family = poisson,
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

test_that("poisson: zero-inflation fixed effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ mined,
    ziformula = ~ mined,
    family = poisson,
    data = Salamanders,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ mined,
    ziformula = ~ mined,
    family = poisson,
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

test_that("poisson: zero-inflation offset", {
  zi_offset_data <- transform(
    Salamanders,
    zi_offset = rep(c(-0.2, 0.2), length.out = nrow(Salamanders))
  )

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ mined,
    ziformula = ~ mined + offset(zi_offset),
    family = poisson,
    data = zi_offset_data,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ mined,
    ziformula = ~ mined + offset(zi_offset),
    family = poisson,
    data = zi_offset_data,
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

test_that("poisson: zero-inflation random intercept", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ mined,
    ziformula = ~ 1 + (1 | site),
    family = poisson,
    data = Salamanders,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ mined,
    ziformula = ~ 1 + (1 | site),
    family = poisson,
    data = Salamanders,
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

test_that("poisson: conditional and zero-inflation random effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ mined + (1 | site),
    ziformula = ~ mined + (1 | site),
    family = poisson,
    data = Salamanders,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ mined + (1 | site),
    ziformula = ~ mined + (1 | site),
    family = poisson,
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
  expect_equal(
    VarCorr(m_rtmb),
    VarCorr(m_tmb),
    tolerance = tol_varcorr
  )
})

test_that("poisson: random-only zero inflation", {
  thetazi <- log(0.5)
  theta_map <- factor(NA)

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ mined,
    ziformula = ~ 0 + (1 | site),
    family = poisson,
    data = Salamanders,
    start = list(thetazi = thetazi),
    map = list(thetazi = theta_map),
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ mined,
    ziformula = ~ 0 + (1 | site),
    family = poisson,
    data = Salamanders,
    start = list(thetazi = thetazi),
    map = list(thetazi = theta_map),
    se = FALSE
  )
  m_nozi <- glmmTMB(
    count ~ mined,
    ziformula = ~ 0,
    family = poisson,
    data = Salamanders,
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
  expect_gt(
    abs(as.numeric(logLik(m_tmb)) - as.numeric(logLik(m_nozi))),
    1
  )
})

test_that("poisson: random-only ZI simulation generates structural zeros", {
  simulation_data <- transform(Salamanders, count = 20L)
  beta <- log(20)
  thetazi <- log(0.5)
  fixed_map <- factor(NA)

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ 1,
    ziformula = ~ 0 + (1 | site),
    family = poisson,
    data = simulation_data,
    start = list(beta = beta, thetazi = thetazi),
    map = list(beta = fixed_map, thetazi = fixed_map),
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ 1,
    ziformula = ~ 0 + (1 | site),
    family = poisson,
    data = simulation_data,
    start = list(beta = beta, thetazi = thetazi),
    map = list(beta = fixed_map, thetazi = fixed_map),
    se = FALSE
  )

  set.seed(2004)
  sim_rtmb <- m_rtmb$obj$simulate(complete = TRUE)$yobs
  set.seed(2004)
  sim_tmb <- m_tmb$obj$simulate(complete = TRUE)$yobs

  expect_gt(sum(sim_rtmb == 0), nrow(simulation_data) / 4)
  expect_gt(sum(sim_tmb == 0), nrow(simulation_data) / 4)
  expect_true(any(sim_rtmb > 0))
  expect_true(any(sim_tmb > 0))
})

test_that("poisson: weighted zero-inflation model", {
  weighted_data <- transform(
    Salamanders,
    w = rep(c(1, 2), length.out = nrow(Salamanders))
  )

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ mined,
    ziformula = ~ mined,
    family = poisson,
    weights = w,
    data = weighted_data,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ mined,
    ziformula = ~ mined,
    family = poisson,
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
  expect_equal(
    fixef(m_rtmb)$zi,
    fixef(m_tmb)$zi,
    tolerance = tol_fixef
  )
})

test_that("poisson: sparse conditional and ZI matrices", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ mined + spp,
    ziformula = ~ mined,
    family = poisson,
    data = Salamanders,
    sparseX = c(cond = TRUE, zi = TRUE),
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ mined + spp,
    ziformula = ~ mined,
    family = poisson,
    data = Salamanders,
    sparseX = c(cond = TRUE, zi = TRUE),
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

test_that("poisson: simulation under RTMB backend", {
  local_useRTMB(TRUE)
  model <- glmmTMB(
    count ~ mined + (1 | site),
    family = poisson,
    data = Salamanders,
    se = FALSE
  )

  set.seed(2002)
  simulated <- model$obj$simulate(complete = TRUE)$yobs

  expect_length(simulated, nrow(Salamanders))
  expect_true(all(is.finite(simulated)))
  expect_true(all(simulated >= 0))
  expect_true(all(simulated == floor(simulated)))
})

test_that("poisson: zero-inflated simulation generates structural zeros", {
  local_useRTMB(TRUE)
  model <- glmmTMB(
    count ~ 1,
    ziformula = ~ 1,
    family = poisson,
    data = Salamanders,
    start = list(beta = log(20), betazi = qlogis(0.5)),
    map = list(beta = factor(NA), betazi = factor(NA)),
    se = FALSE
  )

  set.seed(2003)
  simulated <- model$obj$simulate(complete = TRUE)$yobs

  expect_length(simulated, nrow(Salamanders))
  expect_true(all(is.finite(simulated)))
  expect_true(any(simulated == 0))
  expect_true(any(simulated > 0))
})

test_that("poisson: homogeneous AR1 covariance", {
  theta <- c(log(0.4), 0.3)
  theta_map <- factor(rep(NA, length(theta)))

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ time + ar1(0 + time_fac | group),
    family = poisson,
    data = poisson_panel,
    start = list(theta = theta),
    map = list(theta = theta_map),
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ time + ar1(0 + time_fac | group),
    family = poisson,
    data = poisson_panel,
    start = list(theta = theta),
    map = list(theta = theta_map),
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

test_that("poisson: heterogeneous AR1 covariance", {
  theta <- c(rep(log(0.4), 4), 0.3)
  theta_map <- factor(rep(NA, length(theta)))

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ time + hetar1(0 + time_fac | group),
    family = poisson,
    data = poisson_panel,
    start = list(theta = theta),
    map = list(theta = theta_map),
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ time + hetar1(0 + time_fac | group),
    family = poisson,
    data = poisson_panel,
    start = list(theta = theta),
    map = list(theta = theta_map),
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

test_that("poisson: Ornstein-Uhlenbeck covariance", {
  theta <- c(log(0.4), log(0.7))
  theta_map <- factor(rep(NA, length(theta)))

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ time + ou(0 + time_num | group),
    family = poisson,
    data = poisson_panel,
    start = list(theta = theta),
    map = list(theta = theta_map),
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ time + ou(0 + time_num | group),
    family = poisson,
    data = poisson_panel,
    start = list(theta = theta),
    map = list(theta = theta_map),
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

test_that("poisson: heterogeneous compound-symmetry covariance", {
  theta <- c(rep(log(0.4), 4), 0)
  theta_map <- factor(rep(NA, length(theta)))

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ time + cs(0 + time_fac | group),
    family = poisson,
    data = poisson_panel,
    start = list(theta = theta),
    map = list(theta = theta_map),
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ time + cs(0 + time_fac | group),
    family = poisson,
    data = poisson_panel,
    start = list(theta = theta),
    map = list(theta = theta_map),
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

test_that("poisson: homogeneous compound-symmetry covariance", {
  theta <- c(log(0.4), 0)
  theta_map <- factor(rep(NA, length(theta)))

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ time + homcs(0 + time_fac | group),
    family = poisson,
    data = poisson_panel,
    start = list(theta = theta),
    map = list(theta = theta_map),
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ time + homcs(0 + time_fac | group),
    family = poisson,
    data = poisson_panel,
    start = list(theta = theta),
    map = list(theta = theta_map),
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

test_that("poisson: heterogeneous Toeplitz covariance", {
  theta <- c(rep(log(0.4), 4), rep(0.2, 3))
  theta_map <- factor(rep(NA, length(theta)))

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ time + toep(0 + time_fac | group),
    family = poisson,
    data = poisson_panel,
    start = list(theta = theta),
    map = list(theta = theta_map),
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ time + toep(0 + time_fac | group),
    family = poisson,
    data = poisson_panel,
    start = list(theta = theta),
    map = list(theta = theta_map),
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

test_that("poisson: homogeneous Toeplitz covariance", {
  theta <- c(log(0.4), rep(0.2, 3))
  theta_map <- factor(rep(NA, length(theta)))

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ time + homtoep(0 + time_fac | group),
    family = poisson,
    data = poisson_panel,
    start = list(theta = theta),
    map = list(theta = theta_map),
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ time + homtoep(0 + time_fac | group),
    family = poisson,
    data = poisson_panel,
    start = list(theta = theta),
    map = list(theta = theta_map),
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

test_that("poisson: proportional covariance", {
  matrix_names <- c("(Intercept)", "time")
  proportional_matrix <- diag(2)
  dimnames(proportional_matrix) <- list(matrix_names, matrix_names)

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ time + propto(time | group, proportional_matrix),
    family = poisson,
    data = poisson_panel,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ time + propto(time | group, proportional_matrix),
    family = poisson,
    data = poisson_panel,
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

test_that("poisson: fixed equal-to covariance", {
  matrix_names <- c("(Intercept)", "time")
  fixed_covariance <- matrix(c(0.25, 0.02, 0.02, 0.04), 2, 2)
  dimnames(fixed_covariance) <- list(matrix_names, matrix_names)

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ time + equalto(time | group, fixed_covariance),
    family = poisson,
    data = poisson_panel,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ time + equalto(time | group, fixed_covariance),
    family = poisson,
    data = poisson_panel,
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

test_that("poisson: predict with standard errors", {
  newdata <- Salamanders[seq_len(10), ]
  aggregate <- Salamanders$mined

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ mined + (1 | site),
    family = poisson,
    data = Salamanders,
    se = FALSE
  )
  pred_rtmb <- predict(m_rtmb, se.fit = TRUE)
  link_rtmb <- predict(m_rtmb, type = "link", se.fit = TRUE)
  new_rtmb <- predict(m_rtmb, newdata = newdata, se.fit = TRUE)
  agg_rtmb <- predict(m_rtmb, aggregate = aggregate, se.fit = TRUE)
  cov_rtmb <- predict(m_rtmb, newdata = newdata, se.fit = TRUE,
                      cov.fit = TRUE)
  lat_rtmb <- predict(m_rtmb, type = "latent", se.fit = TRUE)

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ mined + (1 | site),
    family = poisson,
    data = Salamanders,
    se = FALSE
  )
  pred_tmb <- predict(m_tmb, se.fit = TRUE)
  link_tmb <- predict(m_tmb, type = "link", se.fit = TRUE)
  new_tmb <- predict(m_tmb, newdata = newdata, se.fit = TRUE)
  agg_tmb <- predict(m_tmb, aggregate = aggregate, se.fit = TRUE)
  cov_tmb <- predict(m_tmb, newdata = newdata, se.fit = TRUE,
                     cov.fit = TRUE)
  lat_tmb <- predict(m_tmb, type = "latent", se.fit = TRUE)

  expect_equal(pred_rtmb$fit, pred_tmb$fit, tolerance = tol_fixef)
  expect_equal(pred_rtmb$se.fit, pred_tmb$se.fit, tolerance = tol_fixef)
  expect_equal(link_rtmb$fit, link_tmb$fit, tolerance = tol_fixef)
  expect_equal(link_rtmb$se.fit, link_tmb$se.fit, tolerance = tol_fixef)
  expect_equal(new_rtmb$fit, new_tmb$fit, tolerance = tol_fixef)
  expect_equal(new_rtmb$se.fit, new_tmb$se.fit, tolerance = tol_fixef)
  expect_equal(agg_rtmb$fit, agg_tmb$fit, tolerance = tol_fixef)
  expect_equal(agg_rtmb$se.fit, agg_tmb$se.fit, tolerance = tol_fixef)
  expect_equal(cov_rtmb$fit, cov_tmb$fit, tolerance = tol_fixef)
  expect_equal(cov_rtmb$se.fit, cov_tmb$se.fit, tolerance = tol_fixef)
  expect_equal(cov_rtmb$cov.fit, cov_tmb$cov.fit, tolerance = tol_fixef)
  expect_equal(lat_rtmb$fit, lat_tmb$fit, tolerance = tol_fixef)
  expect_equal(lat_rtmb$se.fit, lat_tmb$se.fit, tolerance = tol_fixef)
})

test_that("poisson: zero-inflated predict with standard errors", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ mined + (1 | site),
    ziformula = ~ mined,
    family = poisson,
    data = Salamanders,
    se = FALSE
  )
  resp_rtmb <- predict(m_rtmb, se.fit = TRUE)
  cond_rtmb <- predict(m_rtmb, type = "conditional", se.fit = TRUE)
  zprob_rtmb <- predict(m_rtmb, type = "zprob", se.fit = TRUE)
  zlink_rtmb <- predict(m_rtmb, type = "zlink", se.fit = TRUE)

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ mined + (1 | site),
    ziformula = ~ mined,
    family = poisson,
    data = Salamanders,
    se = FALSE
  )
  resp_tmb <- predict(m_tmb, se.fit = TRUE)
  cond_tmb <- predict(m_tmb, type = "conditional", se.fit = TRUE)
  zprob_tmb <- predict(m_tmb, type = "zprob", se.fit = TRUE)
  zlink_tmb <- predict(m_tmb, type = "zlink", se.fit = TRUE)

  expect_equal(resp_rtmb$fit, resp_tmb$fit, tolerance = tol_fixef)
  expect_equal(resp_rtmb$se.fit, resp_tmb$se.fit, tolerance = tol_fixef)
  expect_equal(cond_rtmb$fit, cond_tmb$fit, tolerance = tol_fixef)
  expect_equal(cond_rtmb$se.fit, cond_tmb$se.fit, tolerance = tol_fixef)
  expect_equal(zprob_rtmb$fit, zprob_tmb$fit, tolerance = tol_fixef)
  expect_equal(zprob_rtmb$se.fit, zprob_tmb$se.fit, tolerance = tol_fixef)
  expect_equal(zlink_rtmb$fit, zlink_tmb$fit, tolerance = tol_fixef)
  expect_equal(zlink_rtmb$se.fit, zlink_tmb$se.fit, tolerance = tol_fixef)
})

test_that("poisson: fixed-effect and random-effect priors", {
  priors <- data.frame(
    prior = c("normal(0, 0.1)", "gamma(1, 2.5)"),
    class = c("fixef", "ranef"),
    coef = c("minedno", "1|site")
  )

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ mined + (1 | site),
    family = poisson,
    data = Salamanders,
    priors = priors,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ mined + (1 | site),
    family = poisson,
    data = Salamanders,
    priors = priors,
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
    getME(m_rtmb, "theta"),
    getME(m_tmb, "theta"),
    tolerance = tol_varcorr
  )
})

test_that("poisson: zero-inflation fixed-effect priors", {
  prior <- data.frame(
    prior = "normal(0, 0.1)",
    class = "fixef_zi",
    coef = "minedno"
  )

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ mined + (1 | site),
    ziformula = ~ mined,
    family = poisson,
    data = Salamanders,
    priors = prior,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ mined + (1 | site),
    ziformula = ~ mined,
    family = poisson,
    data = Salamanders,
    priors = prior,
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
