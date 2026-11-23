context("RTMB one-step-ahead residual support")

skip_if_not_installed("RTMB")

test_that("gaussian OSA residuals match TMB", {
  data("sleepstudy", package = "lme4")

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ Days, data = sleepstudy, family = gaussian,
                    se = FALSE)
  osa_rtmb <- RTMB::oneStepPredict(
    m_rtmb$obj,
    observation.name = "yobs",
    data.term.indicator = "_RTMB_keep_",
    method = "oneStepGaussian",
    discrete = FALSE,
    trace = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ Days, data = sleepstudy, family = gaussian,
                   se = FALSE)
  osa_tmb <- TMB::oneStepPredict(
    m_tmb$obj,
    observation.name = "yobs",
    data.term.indicator = "keep",
    method = "oneStepGaussian",
    discrete = FALSE,
    trace = FALSE
  )

  expect_equal(osa_rtmb$residual, osa_tmb$residual, tolerance = 1e-8)
  expect_equal(osa_rtmb$mean, osa_tmb$mean, tolerance = 1e-8)
  expect_equal(osa_rtmb$sd, osa_tmb$sd, tolerance = 1e-8)
})

test_that("poisson OSA residuals match TMB", {
  data("Salamanders", package = "glmmTMB")
  d <- Salamanders[seq_len(60), ]
  discrete_support <- 0:(max(d$count) + 3L)

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(count ~ mined, data = d, family = poisson, se = FALSE)
  osa_rtmb <- RTMB::oneStepPredict(
    m_rtmb$obj,
    observation.name = "yobs",
    data.term.indicator = "_RTMB_keep_",
    method = "oneStepGeneric",
    discrete = TRUE,
    discreteSupport = discrete_support,
    trace = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(count ~ mined, data = d, family = poisson, se = FALSE)
  osa_tmb <- TMB::oneStepPredict(
    m_tmb$obj,
    observation.name = "yobs",
    data.term.indicator = "keep",
    method = "oneStepGeneric",
    discrete = TRUE,
    discreteSupport = discrete_support,
    trace = FALSE
  )

  expect_equal(osa_rtmb$nll, osa_tmb$nll, tolerance = 1e-8)
  expect_equal(osa_rtmb$Fx, osa_tmb$Fx, tolerance = 1e-8)
  expect_equal(osa_rtmb$px, osa_tmb$px, tolerance = 1e-8)
  expect_equal(osa_rtmb$residual, osa_tmb$residual, tolerance = 1e-8)
})
