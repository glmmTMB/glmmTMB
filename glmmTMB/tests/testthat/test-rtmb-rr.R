## Tests for the RTMB reduced-rank covariance structure

context("RTMB reduced-rank covariance structure")

skip_if_not_installed("RTMB")

tol_logLik <- 1e-5
tol_fixef <- 1e-5
tol_varcorr <- 1e-4

data("spider_long", package = "glmmTMB")
spp_totals <- sort(
  tapply(spider_long$abund, spider_long$Species, sum),
  decreasing = TRUE
)
rr_spider_data <- subset(
  spider_long,
  Species %in% names(spp_totals)[1:4]
)

set.seed(5001)
rr_gaussian_data <- data.frame(
  group = factor(rep(seq_len(30), each = 10)),
  level = factor(rep(seq_len(4), length.out = 300)),
  x = rnorm(300)
)
rr_gaussian_data$response <- 2 + 0.5 * rr_gaussian_data$x +
  rnorm(nrow(rr_gaussian_data))

test_that("Poisson rank-one rr covariance matches TMB", {
  formula <- abund ~ Species + rr(0 + Species | id, d = 1)

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    formula,
    family = poisson,
    data = rr_spider_data,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    formula,
    family = poisson,
    data = rr_spider_data,
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
  expect_equal(
    unname(m_rtmb$obj$report()$fact_load[[1]]),
    unname(m_tmb$obj$report()$fact_load[[1]]),
    tolerance = tol_varcorr
  )
})

test_that("Gaussian rank-two rr covariance matches TMB", {
  formula <- response ~ x + rr(0 + level | group, d = 2)

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    formula,
    family = gaussian,
    data = rr_gaussian_data,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    formula,
    family = gaussian,
    data = rr_gaussian_data,
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
  expect_equal(
    unname(m_rtmb$obj$report()$fact_load[[1]]),
    unname(m_tmb$obj$report()$fact_load[[1]]),
    tolerance = tol_varcorr
  )

  eigenvalues <- eigen(VarCorr(m_rtmb)$cond$group)$values
  expect_equal(
    utils::tail(eigenvalues, 2L),
    c(0, 0),
    tolerance = tol_varcorr
  )
})

test_that("rr covariance works in the zero-inflation model", {
  thetazi <- rep(0.15, nlevels(Salamanders$spp))
  theta_map <- factor(rep(NA, length(thetazi)))

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    count ~ mined,
    ziformula = ~ 1 + rr(0 + spp | site, d = 1),
    family = poisson,
    data = Salamanders,
    start = list(thetazi = thetazi),
    map = list(thetazi = theta_map),
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    count ~ mined,
    ziformula = ~ 1 + rr(0 + spp | site, d = 1),
    family = poisson,
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
    fixef(m_rtmb),
    fixef(m_tmb),
    tolerance = tol_fixef
  )
  expect_equal(
    VarCorr(m_rtmb),
    VarCorr(m_tmb),
    tolerance = tol_varcorr
  )
})

test_that("rr covariance works in the dispersion model", {
  data("sleepstudy", package = "lme4")
  sleepstudy$DaysFac <- factor(sleepstudy$Days)
  thetadisp <- rep(0.08, nlevels(sleepstudy$DaysFac))
  theta_map <- factor(rep(NA, length(thetadisp)))

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    Reaction ~ Days,
    dispformula = ~ 1 + rr(0 + DaysFac | Subject, d = 1),
    family = gaussian,
    data = sleepstudy,
    start = list(thetadisp = thetadisp),
    map = list(thetadisp = theta_map),
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    Reaction ~ Days,
    dispformula = ~ 1 + rr(0 + DaysFac | Subject, d = 1),
    family = gaussian,
    data = sleepstudy,
    start = list(thetadisp = thetadisp),
    map = list(thetadisp = theta_map),
    se = FALSE
  )

  expect_equal(
    as.numeric(logLik(m_rtmb)),
    as.numeric(logLik(m_tmb)),
    tolerance = tol_logLik
  )
  expect_equal(
    fixef(m_rtmb),
    fixef(m_tmb),
    tolerance = tol_fixef
  )
  expect_equal(
    VarCorr(m_rtmb),
    VarCorr(m_tmb),
    tolerance = tol_varcorr
  )
})

test_that("rr simulation works under RTMB", {
  theta <- rep(0.2, 7L)
  theta_map <- factor(rep(NA, length(theta)))

  local_useRTMB(TRUE)
  model <- glmmTMB(
    response ~ x + rr(0 + level | group, d = 2),
    family = gaussian,
    data = rr_gaussian_data,
    start = list(theta = theta),
    map = list(theta = theta_map),
    se = FALSE
  )

  set.seed(5002)
  simulation <- model$obj$simulate(complete = TRUE)$yobs

  expect_length(simulation, nrow(rr_gaussian_data))
  expect_true(all(is.finite(simulation)))
})
