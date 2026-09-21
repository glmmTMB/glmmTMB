## Test Cases for RTMB Gaussian Family
## Tests follow this pattern: fit with both the RTMB and legacy TMB
## backends and compare results (logLik, fixed effects, etc.)

## Example:
# remotes::install_github("glmmTMB/glmmTMB/glmmTMB", ref = "portRTMB")

#  library(glmmTMB)
# local_useRTMB(TRUE)
# m1 <- glmmTMB(count ~ mined + (1|site),
#                family=poisson, data=Salamanders)
# logLik(m1)
# sim <- m1$obj$simulate(complete=TRUE)

# glmmTMB::useRTMB(FALSE)
# m1 <- glmmTMB(count ~ mined + (1|site),
#                family=poisson, data=Salamanders)
# logLik(m1)



##SLEEP STUDY DATASET

context("RTMB Gaussian backend")

skip_if_not_installed("RTMB")

data("sleepstudy", package = "lme4")

tol_logLik <- 1e-5
tol_fixef <- 1e-5
tol_varcorr <- 1e-4

test_that("gaussian: fixed conditional effects only", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ Days, family = gaussian, data = sleepstudy,
                     se = FALSE)

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ Days, family = gaussian, data = sleepstudy,
                    se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(unname(fixef(m_rtmb)$cond), unname(fixef(m_tmb)$cond),
               tolerance = tol_fixef)
})

test_that("gaussian: fixed effects with offset", {
  sleepstudy$off <- 0.1 * sleepstudy$Days

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ Days + offset(off), family = gaussian,
                     data = sleepstudy, se = FALSE)

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ Days + offset(off), family = gaussian,
                    data = sleepstudy, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
})

test_that("gaussian: weighted observations (runif weights)", {
  set.seed(102)
  sleepstudy$w <- runif(nrow(sleepstudy), 0.5, 2)

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ Days, family = gaussian, data = sleepstudy,
                     weights = w, se = FALSE)

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ Days, family = gaussian, data = sleepstudy,
                    weights = w, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
})

test_that("gaussian: weighted fixed effects (alternating weights)", {
  sleepstudy$w <- rep(c(1, 2), length.out = nrow(sleepstudy))

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ Days, family = gaussian, data = sleepstudy,
                     weights = w, se = FALSE)

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ Days, family = gaussian, data = sleepstudy,
                    weights = w, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = tol_fixef)
})

test_that("gaussian: fixed dispersion formula", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ Days, dispformula = ~ Days,
                     family = gaussian, data = sleepstudy, se = FALSE)

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ Days, dispformula = ~ Days,
                    family = gaussian, data = sleepstudy, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = tol_fixef)
  expect_equal(unname(fixef(m_rtmb)$disp), unname(fixef(m_tmb)$disp),
               tolerance = tol_fixef)
})

test_that("gaussian: random effects in dispersion formula", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ Days, dispformula = ~ 1 + (1 | Subject),
                     family = gaussian, data = sleepstudy, se = FALSE)

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ Days, dispformula = ~ 1 + (1 | Subject),
                    family = gaussian, data = sleepstudy, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = tol_fixef)
  expect_equal(VarCorr(m_rtmb), VarCorr(m_tmb), tolerance = tol_varcorr)
})

test_that("gaussian: fixed ZI formula (hurdle, introduced zeros)", {
  set.seed(101)
  sleepstudy$Reaction[sample(nrow(sleepstudy), 5)] <- 0

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ Days, ziformula = ~ Days,
                     family = gaussian, data = sleepstudy, se = FALSE)

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ Days, ziformula = ~ Days,
                    family = gaussian, data = sleepstudy, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(unname(fixef(m_rtmb)$zi), unname(fixef(m_tmb)$zi),
               tolerance = tol_fixef)
})

test_that("gaussian: zero-inflation fixed effects match TMB backend (no induced zeros)", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ Days, ziformula = ~ Days,
                     family = gaussian, data = sleepstudy, se = FALSE)

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ Days, ziformula = ~ Days,
                    family = gaussian, data = sleepstudy, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = tol_fixef)
  expect_equal(fixef(m_rtmb)$zi, fixef(m_tmb)$zi, tolerance = tol_fixef)
})

test_that("gaussian: random-only ZI formula (~0 + RE, induced zeros)", {
  set.seed(104)
  sleepstudy$Reaction[sample(nrow(sleepstudy), nrow(sleepstudy) / 2)] <- 0
  thetazi <- log(0.5)
  theta_map <- factor(NA)

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ Days, ziformula = ~ 0 + (1 | Subject),
                     family = gaussian, data = sleepstudy,
                     start = list(thetazi = thetazi),
                     map = list(thetazi = theta_map), se = FALSE)

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ Days, ziformula = ~ 0 + (1 | Subject),
                    family = gaussian, data = sleepstudy,
                    start = list(thetazi = thetazi),
                    map = list(thetazi = theta_map), se = FALSE)
  m_nozi <- glmmTMB(Reaction ~ Days, ziformula = ~ 0,
                     family = gaussian, data = sleepstudy, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(VarCorr(m_rtmb), VarCorr(m_tmb), tolerance = tol_varcorr)
  expect_gt(
    abs(as.numeric(logLik(m_tmb)) - as.numeric(logLik(m_nozi))),
    1
  )
})

test_that("gaussian: ZI intercept + random effects match TMB backend (no induced zeros)", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ Days, ziformula = ~ 1 + (1 | Subject),
                     family = gaussian, data = sleepstudy, se = FALSE)

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ Days, ziformula = ~ 1 + (1 | Subject),
                    family = gaussian, data = sleepstudy, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = tol_fixef)
  expect_equal(VarCorr(m_rtmb), VarCorr(m_tmb), tolerance = tol_varcorr)
})

test_that("gaussian: dZI simulation generates structural zeros", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    Reaction ~ Days,
    ziformula = ~ 1,
    family = gaussian,
    data = sleepstudy,
    start = list(betazi = qlogis(0.5)),
    map = list(betazi = factor(NA)),
    se = FALSE
  )

  set.seed(1001)
  sim <- m_rtmb$obj$simulate(complete = TRUE)$yobs

  expect_length(sim, nrow(sleepstudy))
  expect_true(all(is.finite(sim)))
  expect_true(any(sim == 0))
  expect_true(any(sim != 0))
})

test_that("gaussian: homogeneous AR1 covariance", {
  ar1_dat <- transform(sleepstudy, DaysFac = factor(Days))

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    Reaction ~ Days + ar1(0 + DaysFac | Subject),
    family = gaussian,
    data = ar1_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    Reaction ~ Days + ar1(0 + DaysFac | Subject),
    family = gaussian,
    data = ar1_dat,
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

test_that("gaussian: AR1 covariance matrix follows phi distance", {
  ar1_dat <- transform(sleepstudy, DaysFac = factor(Days))

  local_useRTMB(TRUE)
  model <- glmmTMB(
    Reaction ~ Days + ar1(0 + DaysFac | Subject),
    family = gaussian,
    data = ar1_dat,
    se = FALSE
  )

  vc <- VarCorr(model)$cond$Subject
  correlation <- attr(vc, "correlation")
  phi <- correlation[1L, 2L]
  expected <- phi^abs(
    outer(seq_len(nrow(correlation)), seq_len(ncol(correlation)), "-")
  )
  dimnames(expected) <- dimnames(correlation)

  expect_equal(correlation, expected, tolerance = tol_varcorr)
})

test_that("gaussian: AR1 works in zero-inflation model", {
  ar1_dat <- transform(sleepstudy, DaysFac = factor(Days))
  set.seed(301)
  ar1_dat$Reaction[sample(nrow(ar1_dat), 20)] <- 0
  rho <- 0.4
  thetazi <- c(log(0.5), rho / sqrt(1 - rho^2))
  theta_map <- factor(c(NA, NA))

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    Reaction ~ Days,
    ziformula = ~ 1 + ar1(0 + DaysFac | Subject),
    family = gaussian,
    data = ar1_dat,
    start = list(thetazi = thetazi),
    map = list(thetazi = theta_map),
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    Reaction ~ Days,
    ziformula = ~ 1 + ar1(0 + DaysFac | Subject),
    family = gaussian,
    data = ar1_dat,
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

test_that("gaussian: heterogeneous AR1 covariance", {
  set.seed(501)

  ng <- 60L
  nt <- 3L
  rho <- 0.45
  component_sd <- c(0.7, 1.1, 1.6)
  correlation <- rho^abs(
    outer(seq_len(nt), seq_len(nt), "-")
  )
  Sigma <- outer(component_sd, component_sd) * correlation
  random_effects <- MASS::mvrnorm(ng, rep(0, nt), Sigma)

  hetar1_dat <- data.frame(
    group = factor(rep(seq_len(ng), each = nt)),
    time = factor(rep(seq_len(nt), times = ng)),
    time_num = rep(seq_len(nt) - 1, times = ng)
  )
  hetar1_dat$y <- 2 + 0.35 * hetar1_dat$time_num +
    as.vector(t(random_effects)) +
    rnorm(nrow(hetar1_dat), sd = 0.4)

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ time_num + hetar1(0 + time | group),
    family = gaussian,
    data = hetar1_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ time_num + hetar1(0 + time | group),
    family = gaussian,
    data = hetar1_dat,
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

test_that("gaussian: hetar1 reports component-specific SDs", {
  ar1_dat <- transform(sleepstudy, DaysFac = factor(Days))
  component_sd <- seq(0.5, 1.4, length.out = nlevels(ar1_dat$DaysFac))
  rho <- 0.4
  theta <- c(
    log(component_sd),
    rho / sqrt(1 - rho^2)
  )

  local_useRTMB(TRUE)
  model <- glmmTMB(
    Reaction ~ Days + hetar1(0 + DaysFac | Subject),
    family = gaussian,
    data = ar1_dat,
    start = list(theta = theta),
    map = list(theta = factor(rep(NA, length(theta)))),
    se = FALSE
  )

  vc <- VarCorr(model)$cond$Subject
  correlation <- attr(vc, "correlation")
  reported_sd <- attr(vc, "stddev")
  expected_corr <- rho^abs(
    outer(
      seq_len(nrow(correlation)),
      seq_len(ncol(correlation)),
      "-"
    )
  )
  dimnames(expected_corr) <- dimnames(correlation)

  expect_equal(unname(reported_sd), component_sd, tolerance = tol_varcorr)
  expect_equal(correlation, expected_corr, tolerance = tol_varcorr)
})

test_that("gaussian: hetar1 works in zero-inflation model", {
  ar1_dat <- transform(sleepstudy, DaysFac = factor(Days))
  set.seed(502)
  ar1_dat$Reaction[sample(nrow(ar1_dat), 20)] <- 0

  component_sd <- rep(0.5, nlevels(ar1_dat$DaysFac))
  rho <- 0.4
  thetazi <- c(
    log(component_sd),
    rho / sqrt(1 - rho^2)
  )
  theta_map <- factor(rep(NA, length(thetazi)))

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    Reaction ~ Days,
    ziformula = ~ 1 + hetar1(0 + DaysFac | Subject),
    family = gaussian,
    data = ar1_dat,
    start = list(thetazi = thetazi),
    map = list(thetazi = theta_map),
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    Reaction ~ Days,
    ziformula = ~ 1 + hetar1(0 + DaysFac | Subject),
    family = gaussian,
    data = ar1_dat,
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

test_that("gaussian: Ornstein-Uhlenbeck covariance", {
  set.seed(601)

  ng <- 40L
  times <- c(0, 0.5, 2, 5)
  nt <- length(times)
  decay <- 0.7
  process_sd <- 1.2
  correlation <- exp(
    -decay * abs(outer(times, times, "-"))
  )
  Sigma <- process_sd^2 * correlation
  random_effects <- MASS::mvrnorm(ng, rep(0, nt), Sigma)

  ou_dat <- data.frame(
    group = factor(rep(seq_len(ng), each = nt)),
    time = glmmTMB::numFactor(rep(times, times = ng)),
    time_num = rep(times, times = ng)
  )
  ou_dat$y <- 1 + 0.2 * ou_dat$time_num +
    as.vector(t(random_effects))

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y ~ time_num + ou(0 + time | group),
    family = gaussian,
    dispformula = ~ 0,
    data = ou_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y ~ time_num + ou(0 + time | group),
    family = gaussian,
    dispformula = ~ 0,
    data = ou_dat,
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

test_that("gaussian: OU matrix follows continuous time distances", {
  set.seed(602)

  times <- c(0, 0.5, 2, 5)
  ng <- 20L
  process_sd <- 1.2
  decay <- 0.7
  theta <- c(log(process_sd), log(decay))
  ou_dat <- data.frame(
    group = factor(rep(seq_len(ng), each = length(times))),
    time = glmmTMB::numFactor(rep(times, times = ng)),
    y = rnorm(ng * length(times))
  )

  local_useRTMB(TRUE)
  model <- glmmTMB(
    y ~ ou(0 + time | group),
    family = gaussian,
    data = ou_dat,
    start = list(theta = theta),
    map = list(theta = factor(c(NA, NA))),
    se = FALSE
  )

  vc <- VarCorr(model)$cond$group
  reported_sd <- attr(vc, "stddev")
  correlation <- attr(vc, "correlation")
  expected_corr <- exp(
    -decay * abs(outer(times, times, "-"))
  )
  dimnames(expected_corr) <- dimnames(correlation)

  expect_equal(
    unname(reported_sd),
    rep(process_sd, length(times)),
    tolerance = tol_varcorr
  )
  expect_equal(
    correlation,
    expected_corr,
    tolerance = tol_varcorr
  )
})

test_that("gaussian: OU works in zero-inflation model", {
  ou_dat <- transform(
    sleepstudy,
    DaysNum = glmmTMB::numFactor(Days)
  )
  set.seed(603)
  ou_dat$Reaction[sample(nrow(ou_dat), 20)] <- 0

  process_sd <- 0.5
  decay <- 0.7
  thetazi <- c(log(process_sd), log(decay))
  theta_map <- factor(c(NA, NA))

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    Reaction ~ Days,
    ziformula = ~ 1 + ou(0 + DaysNum | Subject),
    family = gaussian,
    data = ou_dat,
    start = list(thetazi = thetazi),
    map = list(thetazi = theta_map),
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    Reaction ~ Days,
    ziformula = ~ 1 + ou(0 + DaysNum | Subject),
    family = gaussian,
    data = ou_dat,
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

test_that("gaussian: single random intercept (cond RE)", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ Days + (1 | Subject), family = gaussian,
                     data = sleepstudy, se = FALSE)

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ Days + (1 | Subject), family = gaussian,
                    data = sleepstudy, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(as.numeric(VarCorr(m_rtmb)$cond$Subject),
               as.numeric(VarCorr(m_tmb)$cond$Subject),
               tolerance = tol_varcorr)
})

test_that("gaussian: correlated random slope (us covstruct)", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ Days + (Days | Subject), family = gaussian,
                     data = sleepstudy, se = FALSE)

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ Days + (Days | Subject), family = gaussian,
                    data = sleepstudy, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
})

test_that("gaussian: multiple random-effect terms", {
  set.seed(103)
  sleepstudy$grp2 <- factor(sample(1:5, nrow(sleepstudy), replace = TRUE))

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ Days + (1 | Subject) + (1 | grp2),
                     family = gaussian, data = sleepstudy, se = FALSE)

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ Days + (1 | Subject) + (1 | grp2),
                    family = gaussian, data = sleepstudy, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
})

test_that("gaussian: diag covstruct random effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ Days + diag(Days | Subject),
                     family = gaussian, data = sleepstudy, se = FALSE)

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ Days + diag(Days | Subject),
                    family = gaussian, data = sleepstudy, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
})

test_that("gaussian: simulate() works under RTMB backend", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ Days + (1 | Subject), family = gaussian,
                     data = sleepstudy, se = FALSE)
  sim <- m_rtmb$obj$simulate(complete = TRUE)

  expect_true(is.list(sim))
  expect_equal(length(sim$yobs), nrow(sleepstudy))
  expect_true(all(is.finite(sim$yobs)))
})






## Salamander Dataset

data("Salamanders", package = "glmmTMB")
Salamanders$lcount <- log1p(Salamanders$count)  ## continuous response for gaussian

test_that("gaussian (Salamanders): fixed effects with binary factor predictor", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(lcount ~ mined, family = gaussian, data = Salamanders,
                     se = FALSE)

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(lcount ~ mined, family = gaussian, data = Salamanders,
                    se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(unname(fixef(m_rtmb)$cond), unname(fixef(m_tmb)$cond),
               tolerance = tol_fixef)
})

# test_that("gaussian (Salamanders): fixed effects with multi-level factor predictor", {
#   local_useRTMB(TRUE)
#   m_rtmb <- glmmTMB(lcount ~ mined + spp, family = gaussian,
#                      data = Salamanders, se = FALSE)
#
#   glmmTMB::useRTMB(FALSE)
#   m_tmb <- glmmTMB(lcount ~ mined + spp, family = gaussian,
#                     data = Salamanders, se = FALSE)
#
#   expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
#                tolerance = tol_logLik)
#   expect_equal(unname(fixef(m_rtmb)$cond), unname(fixef(m_tmb)$cond),
#                tolerance = tol_fixef)
# })

test_that("gaussian (Salamanders): single random intercept by site", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(lcount ~ mined + (1 | site), family = gaussian,
                     data = Salamanders, se = FALSE)

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(lcount ~ mined + (1 | site), family = gaussian,
                    data = Salamanders, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(as.numeric(VarCorr(m_rtmb)$cond$site),
               as.numeric(VarCorr(m_tmb)$cond$site),
               tolerance = tol_varcorr)
})

test_that("gaussian (Salamanders): crossed random intercepts (site + spp)", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(lcount ~ mined + (1 | site) + (1 | spp),
                     family = gaussian, data = Salamanders, se = FALSE)

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(lcount ~ mined + (1 | site) + (1 | spp),
                    family = gaussian, data = Salamanders, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
})

test_that("gaussian (Salamanders): nested random slope by site", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(lcount ~ mined + (mined | site), family = gaussian,
                     data = Salamanders, se = FALSE)

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(lcount ~ mined + (mined | site), family = gaussian,
                    data = Salamanders, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
})

test_that("gaussian (Salamanders): dispersion varying by mined status", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(lcount ~ mined, dispformula = ~ mined,
                     family = gaussian, data = Salamanders, se = FALSE)

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(lcount ~ mined, dispformula = ~ mined,
                    family = gaussian, data = Salamanders, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(unname(fixef(m_rtmb)$disp), unname(fixef(m_tmb)$disp),
               tolerance = tol_fixef)
})

test_that("gaussian (Salamanders): zero-inflation with factor predictor (hurdle)", {
  ## lcount already has structural zeros from count == 0 observations,
  ## so no artificial zero injection needed here
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(lcount ~ mined, ziformula = ~ mined,
                     family = gaussian, data = Salamanders, se = FALSE)

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(lcount ~ mined, ziformula = ~ mined,
                    family = gaussian, data = Salamanders, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(unname(fixef(m_rtmb)$zi), unname(fixef(m_tmb)$zi),
               tolerance = tol_fixef)
})

test_that("gaussian (Salamanders): zero-inflation with random intercept by site", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(lcount ~ mined, ziformula = ~ 1 + (1 | site),
                     family = gaussian, data = Salamanders, se = FALSE)

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(lcount ~ mined, ziformula = ~ 1 + (1 | site),
                    family = gaussian, data = Salamanders, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = tol_fixef)
})

test_that("gaussian (Salamanders): weighted observations", {
  set.seed(201)
  Salamanders$w <- runif(nrow(Salamanders), 0.5, 2)

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(lcount ~ mined, family = gaussian, data = Salamanders,
                     weights = w, se = FALSE)

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(lcount ~ mined, family = gaussian, data = Salamanders,
                    weights = w, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
})

test_that("gaussian (Salamanders): simulate() works under RTMB backend", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(lcount ~ mined + (1 | site), family = gaussian,
                     data = Salamanders, se = FALSE)
  sim <- m_rtmb$obj$simulate(complete = TRUE)

  expect_true(is.list(sim))
  expect_equal(length(sim$yobs), nrow(Salamanders))
  expect_true(all(is.finite(sim$yobs)))
})

test_that("gaussian: cloglog preserves small probabilities and derivatives", {
  local_useRTMB(TRUE)
  eta <- c(-60, -40, -20, 0, 2)
  mu <- -expm1(-exp(eta))
  dat <- data.frame(y = mu + c(0.1, -0.1, 0.1, -0.1, 0.1), eta = eta)
  fit <- glmmTMB(
    y ~ 1 + offset(eta), data = dat, family = gaussian(link = "cloglog"),
    start = list(beta = 0), map = list(beta = factor(NA)), se = FALSE
  )

  ## Compare relative errors so rounding a tiny probability to zero fails.
  expect_equal(as.numeric(predict(fit, type = "response")) / mu,
               rep(1, length(mu)), tolerance = 1e-12)

  obj <- RTMB::MakeADFun(
    function(p) sum(log_inverse_linkfun_rtmb(p$eta, "cloglog")),
    list(eta = eta), silent = TRUE
  )
  expect_equal(obj$fn(), sum(log(mu)), tolerance = 1e-12)
  expected_gradient <- exp(eta - exp(eta)) / mu
  expect_equal(as.vector(obj$gr()) / expected_gradient,
               rep(1, length(eta)), tolerance = 1e-12)
})

test_that("gaussian: supported inverse links match TMB and manual likelihood", {
  link_parameters <- list(
    identity = c(0.3, 0.4),
    log = c(-0.3, 0.4),
    inverse = c(2, 0.3),
    sqrt = c(0.8, 0.2),
    logit = c(-0.4, 0.8),
    probit = c(-0.2, 0.6),
    cloglog = c(-0.4, 0.5)
  )
  x <- seq(-0.5, 0.5, length.out = 80)
  residual_sd <- 0.2

  inverse_link <- function(link, eta) {
    switch(
      link,
      identity = eta,
      log = exp(eta),
      inverse = 1 / eta,
      sqrt = eta^2,
      logit = 1 / (1 + exp(-eta)),
      probit = stats::pnorm(eta),
      cloglog = 1 - exp(-exp(eta))
    )
  }

  for (link in names(link_parameters)) {
    beta <- link_parameters[[link]]
    eta <- beta[1L] + beta[2L] * x
    mu <- inverse_link(link, eta)
    link_data <- data.frame(
      y = mu + 0.05 * sin(seq_along(x)),
      x = x
    )
    start <- list(beta = beta, betadisp = log(residual_sd))

    local_useRTMB(TRUE)
    m_rtmb <- glmmTMB(
      y ~ x,
      family = gaussian(link = link),
      data = link_data,
      start = start,
      se = FALSE
    )

    glmmTMB::useRTMB(FALSE)
    m_tmb <- glmmTMB(
      y ~ x,
      family = gaussian(link = link),
      data = link_data,
      start = start,
      se = FALSE
    )

    expect_equal(
      as.numeric(logLik(m_rtmb)),
      as.numeric(logLik(m_tmb)),
      tolerance = tol_logLik,
      info = link
    )
    expect_equal(
      fixef(m_rtmb)$cond,
      fixef(m_tmb)$cond,
      tolerance = tol_fixef,
      info = link
    )

    simulation <- m_rtmb$obj$simulate(complete = TRUE)$yobs
    expect_length(simulation, nrow(link_data))
    expect_true(all(is.finite(simulation)), info = link)
  }
})







## ChickWeight Dataset

data("ChickWeight", package = "datasets")
chick_dat <- transform(ChickWeight,
                       Chick = factor(Chick),
                       Diet = factor(Diet))
chick_dat$off <- 0.05 * chick_dat$Time
chick_dat$w <- rep(c(1, 1.5), length.out = nrow(chick_dat))
chick_dat$grp2 <- factor(as.integer(chick_dat$Chick) %% 6)

test_that("gaussian ChickWeight: fixed conditional effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(weight ~ Time + Diet, family = gaussian,
                    data = chick_dat, se = FALSE)

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(weight ~ Time + Diet, family = gaussian,
                   data = chick_dat, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = tol_fixef)
})

test_that("gaussian ChickWeight: fixed effects with offset", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(weight ~ Time + Diet + offset(off), family = gaussian,
                    data = chick_dat, se = FALSE)

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(weight ~ Time + Diet + offset(off), family = gaussian,
                   data = chick_dat, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
})

test_that("gaussian ChickWeight: weighted observations", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(weight ~ Time + Diet, family = gaussian,
                    data = chick_dat, weights = w, se = FALSE)

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(weight ~ Time + Diet, family = gaussian,
                   data = chick_dat, weights = w, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = tol_fixef)
})

test_that("gaussian ChickWeight: fixed dispersion formula", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(weight ~ Time + Diet, dispformula = ~ Time,
                    family = gaussian, data = chick_dat, se = FALSE)

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(weight ~ Time + Diet, dispformula = ~ Time,
                   family = gaussian, data = chick_dat, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(fixef(m_rtmb)$disp, fixef(m_tmb)$disp, tolerance = tol_fixef)
})

test_that("gaussian ChickWeight: conditional random intercept", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(weight ~ Time + Diet + (1 | Chick),
                    family = gaussian, data = chick_dat, se = FALSE)

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(weight ~ Time + Diet + (1 | Chick),
                   family = gaussian, data = chick_dat, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(VarCorr(m_rtmb), VarCorr(m_tmb), tolerance = tol_varcorr)
})

test_that("gaussian ChickWeight: conditional random slope", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(weight ~ Time + Diet + (Time | Chick),
                    family = gaussian, data = chick_dat, se = FALSE)

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(weight ~ Time + Diet + (Time | Chick),
                   family = gaussian, data = chick_dat, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(VarCorr(m_rtmb), VarCorr(m_tmb), tolerance = tol_varcorr)
})

test_that("gaussian ChickWeight: diag covariance random effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(weight ~ Time + Diet + diag(Time | Chick),
                    family = gaussian, data = chick_dat, se = FALSE)

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(weight ~ Time + Diet + diag(Time | Chick),
                   family = gaussian, data = chick_dat, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
})

test_that("gaussian ChickWeight: multiple random-effect terms", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(weight ~ Time + Diet + (1 | Chick) + (1 | grp2),
                    family = gaussian, data = chick_dat, se = FALSE)

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(weight ~ Time + Diet + (1 | Chick) + (1 | grp2),
                   family = gaussian, data = chick_dat, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
})

test_that("gaussian ChickWeight: fixed ZI formula with induced zeros", {
  chick_zi <- chick_dat
  set.seed(201)
  chick_zi$weight[sample(nrow(chick_zi), 10)] <- 0

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(weight ~ Time + Diet, ziformula = ~ Time,
                    family = gaussian, data = chick_zi, se = FALSE)

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(weight ~ Time + Diet, ziformula = ~ Time,
                   family = gaussian, data = chick_zi, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(fixef(m_rtmb)$zi, fixef(m_tmb)$zi, tolerance = tol_fixef)
})

test_that("gaussian ChickWeight: ZI random effects with induced zeros", {
  chick_zi <- chick_dat
  set.seed(202)
  chick_zi$weight[sample(nrow(chick_zi), 10)] <- 0

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(weight ~ Time + Diet, ziformula = ~ 1 + (1 | Chick),
                    family = gaussian, data = chick_zi, se = FALSE)

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(weight ~ Time + Diet, ziformula = ~ 1 + (1 | Chick),
                   family = gaussian, data = chick_zi, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(VarCorr(m_rtmb), VarCorr(m_tmb), tolerance = tol_varcorr)
})





## Existing Gaussian test coverage ported to RTMB/TMB comparisons

test_that("gaussian existing basics: intercept-only random effect", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ 1 + (1 | Subject), data = sleepstudy,
                    se = FALSE)

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ 1 + (1 | Subject), data = sleepstudy,
                   se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = tol_fixef)
  expect_equal(VarCorr(m_rtmb), VarCorr(m_tmb), tolerance = tol_varcorr)
})

test_that("gaussian existing basics: split random intercept and slope", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ Days + (1 | Subject) + (0 + Days | Subject),
                    data = sleepstudy, se = FALSE)

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ Days + (1 | Subject) + (0 + Days | Subject),
                   data = sleepstudy, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = tol_fixef)
  expect_equal(VarCorr(m_rtmb), VarCorr(m_tmb), tolerance = tol_varcorr)
})

test_that("gaussian existing basics: double-bar random effects", {
  data("sleepstudy", package = "lme4")

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ 1 + (Days || Subject), data = sleepstudy,
                    se = FALSE)

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ 1 + (Days || Subject), data = sleepstudy,
                   se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = tol_fixef)
  expect_equal(VarCorr(m_rtmb), VarCorr(m_tmb), tolerance = tol_varcorr)
})

test_that("gaussian existing basics: bar/double-bar bug model", {
  set.seed(1)
  n <- 100
  xdata <- data.frame(
    rfac1 = as.factor(sample(letters[1:10], n, replace = TRUE)),
    rfac2 = as.factor(sample(letters[1:10], n, replace = TRUE)),
    cov = rnorm(n),
    rv = rpois(n, lambda = 2)
  )

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(rv ~ cov + (1 + cov || rfac1) + (1 | rfac2),
                    family = gaussian, data = xdata, se = FALSE)

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(rv ~ cov + (1 + cov || rfac1) + (1 | rfac2),
                   family = gaussian, data = xdata, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = tol_fixef)
  expect_equal(VarCorr(m_rtmb), VarCorr(m_tmb), tolerance = tol_varcorr)
})

test_that("gaussian existing Anova case: fixed dispersion indicator", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ Days + (1 | Subject),
                    dispformula = ~ I(Days > 5), data = sleepstudy,
                    REML = FALSE, se = FALSE)

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ Days + (1 | Subject),
                   dispformula = ~ I(Days > 5), data = sleepstudy,
                   REML = FALSE, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = tol_fixef)
  expect_equal(fixef(m_rtmb)$disp, fixef(m_tmb)$disp, tolerance = tol_fixef)
  expect_equal(VarCorr(m_rtmb), VarCorr(m_tmb), tolerance = tol_varcorr)
})

test_that("gaussian existing methods: Salamanders dispersion by cover", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(count ~ cover, family = gaussian,
                    dispformula = ~ cover, data = Salamanders, se = FALSE)

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(count ~ cover, family = gaussian,
                   dispformula = ~ cover, data = Salamanders, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = tol_fixef)
  expect_equal(fixef(m_rtmb)$disp, fixef(m_tmb)$disp, tolerance = tol_fixef)
})

test_that("gaussian existing methods: mtcars random dispersion effect", {
  mtcars$cyl <- factor(mtcars$cyl)

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(mpg ~ wt, dispformula = ~ 1 + (1 | cyl),
                    data = mtcars, family = gaussian, se = FALSE)

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(mpg ~ wt, dispformula = ~ 1 + (1 | cyl),
                   data = mtcars, family = gaussian, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = tol_fixef)
  expect_equal(VarCorr(m_rtmb), VarCorr(m_tmb), tolerance = tol_varcorr)
})

test_that("gaussian existing varstruc: homdiag random effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ Days + homdiag(Days | Subject),
                    data = sleepstudy, se = FALSE)

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ Days + homdiag(Days | Subject),
                   data = sleepstudy, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = tol_fixef)
  expect_equal(VarCorr(m_rtmb), VarCorr(m_tmb), tolerance = tol_varcorr)
})

test_that("gaussian existing predict: polynomial fixed effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ poly(Days, 3), data = sleepstudy,
                    se = FALSE)

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ poly(Days, 3), data = sleepstudy,
                   se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = tol_fixef)
})

test_that("gaussian existing predict: polynomial with random intercept", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(Reaction ~ (1 | Subject) + poly(Days, 3),
                    data = sleepstudy, se = FALSE)

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(Reaction ~ (1 | Subject) + poly(Days, 3),
                   data = sleepstudy, se = FALSE)

  expect_equal(as.numeric(logLik(m_rtmb)), as.numeric(logLik(m_tmb)),
               tolerance = tol_logLik)
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = tol_fixef)
  expect_equal(VarCorr(m_rtmb), VarCorr(m_tmb), tolerance = tol_varcorr)
})


## Testing covariance structure options

test_that("gaussian: heterogeneous compound-symmetry covariance", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    Reaction ~ Days + cs(Days | Subject),
    family = gaussian,
    data = sleepstudy,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    Reaction ~ Days + cs(Days | Subject),
    family = gaussian,
    data = sleepstudy,
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

test_that("gaussian: homogeneous compound-symmetry covariance", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    Reaction ~ Days + homcs(Days | Subject),
    family = gaussian,
    data = sleepstudy,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    Reaction ~ Days + homcs(Days | Subject),
    family = gaussian,
    data = sleepstudy,
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

test_that("gaussian: heterogeneous Toeplitz covariance", {
  toep_data <- sleepstudy
  toep_data$Reaction <- ifelse(toep_data$Reaction > 250, 1, 0)
  toep_data$Days <- cut(
    toep_data$Days,
    breaks = c(0, 3, 6, 10),
    right = FALSE
  )

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    Reaction ~ toep(0 + Days | Subject),
    family = gaussian,
    data = toep_data,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    Reaction ~ toep(0 + Days | Subject),
    family = gaussian,
    data = toep_data,
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

test_that("gaussian: homogeneous Toeplitz covariance", {
  toep_data <- sleepstudy
  toep_data$Reaction <- ifelse(toep_data$Reaction > 250, 1, 0)
  toep_data$Days <- cut(
    toep_data$Days,
    breaks = c(0, 3, 6, 10),
    right = FALSE
  )

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    Reaction ~ homtoep(0 + Days | Subject),
    family = gaussian,
    data = toep_data,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    Reaction ~ homtoep(0 + Days | Subject),
    family = gaussian,
    data = toep_data,
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

test_that("gaussian: proportional covariance", {
  matrix_names <- c("(Intercept)", "Days")
  proportional_matrix <- diag(2)
  dimnames(proportional_matrix) <- list(matrix_names, matrix_names)

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    Reaction ~ Days + propto(Days | Subject, proportional_matrix),
    family = gaussian,
    data = sleepstudy,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    Reaction ~ Days + propto(Days | Subject, proportional_matrix),
    family = gaussian,
    data = sleepstudy,
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

test_that("gaussian: fixed equal-to covariance", {
  matrix_names <- c("(Intercept)", "Days")
  fixed_covariance <- matrix(c(900, 5, 5, 25), 2, 2)
  dimnames(fixed_covariance) <- list(matrix_names, matrix_names)

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    Reaction ~ Days + equalto(Days | Subject, fixed_covariance),
    family = gaussian,
    data = sleepstudy,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    Reaction ~ Days + equalto(Days | Subject, fixed_covariance),
    family = gaussian,
    data = sleepstudy,
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

test_that("RTMB uses TMB's unstructured theta ordering in dimension four", {
  L <- matrix(
    c(
      1.0, 0.0, 0.0, 0.0,
      0.2, 1.0, 0.0, 0.0,
      0.3, 0.4, 1.0, 0.0,
      0.5, 0.6, 0.7, 1.0
    ),
    nrow = 4,
    byrow = TRUE
  )
  scale <- diag(c(0.5, 1.0, 1.5, 2.0))
  Sigma <- scale %*% tcrossprod(L) %*% scale
  theta <- glmmTMB:::as.theta.vcov(Sigma)

  equalto_term <- list(
    blockSize = 4,
    blockReps = 1,
    blockNumTheta = 10,
    blockCode = structure(15, names = "equalto"),
    fullCor = 1
  )
  equalto_result <- glmmTMB:::termwise_nll(
    numeric(4), theta, equalto_term
  )

  expect_equal(equalto_result$corr, cov2cor(Sigma), tolerance = 1e-12)
  expect_equal(equalto_result$sd, sqrt(diag(Sigma)), tolerance = 1e-12)

  loglambda <- log(2.5)
  propto_term <- equalto_term
  propto_term$blockNumTheta <- 11
  propto_term$blockCode <- structure(11, names = "propto")
  propto_result <- glmmTMB:::termwise_nll(
    numeric(4), c(theta, loglambda), propto_term
  )
  propto_cov <- propto_result$corr *
    tcrossprod(propto_result$sd)

  expect_equal(propto_cov, exp(loglambda) * Sigma, tolerance = 1e-12)
})

test_that("RTMB covariance reports respect full_cor = FALSE", {
  make_term <- function(name, block_size, block_num_theta) {
    list(
      blockSize = block_size,
      blockReps = 1,
      blockNumTheta = block_num_theta,
      blockCode = structure(
        unname(glmmTMB:::.valid_covstruct[[name]]),
        names = name
      ),
      fullCor = 0
    )
  }

  cases <- list(
    us = list(n = 2, theta = c(0, 0, 0.2)),
    cs = list(n = 2, theta = c(0, 0, 0.2)),
    homcs = list(n = 2, theta = c(0, 0.2)),
    toep = list(n = 3, theta = c(0, 0, 0, 0.2, 0.1)),
    homtoep = list(n = 3, theta = c(0, 0.2, 0.1)),
    propto = list(n = 2, theta = c(0, 0, 0.2, 0))
  )

  for (name in names(cases)) {
    case <- cases[[name]]
    term <- make_term(name, case$n, length(case$theta))
    result <- glmmTMB:::termwise_nll(
      numeric(case$n), case$theta, term
    )

    expect_identical(
      result$corr,
      matrix(NaN, 1, 1),
      info = name
    )
  }

  ## The C++ equalto implementation always reports its fixed correlation
  ## matrix, even when full_cor is false.
  equalto_term <- make_term("equalto", 2, 3)
  equalto_result <- glmmTMB:::termwise_nll(
    numeric(2), c(0, 0, 0.2), equalto_term
  )
  expect_equal(dim(equalto_result$corr), c(2, 2))
})

test_that("RTMB full_cor = FALSE reporting works through MakeADFun", {
  local_useRTMB(TRUE)
  fit <- glmmTMB(
    Reaction ~ Days + (Days | Subject),
    family = gaussian,
    data = sleepstudy,
    se = FALSE,
    control = glmmTMBControl(full_cor = FALSE)
  )

  expect_identical(fit$obj$report()$corr[[1]], matrix(NaN, 1, 1))
})

test_that("gaussian: predict with standard errors", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    Reaction ~ Days + (1 | Subject),
    family = gaussian,
    data = sleepstudy,
    se = FALSE
  )
  pred_rtmb <- predict(m_rtmb, se.fit = TRUE)
  link_rtmb <- predict(m_rtmb, type = "link", se.fit = TRUE)
  disp_rtmb <- predict(m_rtmb, type = "disp", se.fit = TRUE)

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    Reaction ~ Days + (1 | Subject),
    family = gaussian,
    data = sleepstudy,
    se = FALSE
  )
  pred_tmb <- predict(m_tmb, se.fit = TRUE)
  link_tmb <- predict(m_tmb, type = "link", se.fit = TRUE)
  disp_tmb <- predict(m_tmb, type = "disp", se.fit = TRUE)

  expect_equal(pred_rtmb$fit, pred_tmb$fit, tolerance = tol_fixef)
  expect_equal(pred_rtmb$se.fit, pred_tmb$se.fit, tolerance = tol_fixef)
  expect_equal(link_rtmb$fit, link_tmb$fit, tolerance = tol_fixef)
  expect_equal(link_rtmb$se.fit, link_tmb$se.fit, tolerance = tol_fixef)
  expect_equal(disp_rtmb$fit, disp_tmb$fit, tolerance = tol_fixef)
  expect_equal(disp_rtmb$se.fit, disp_tmb$se.fit, tolerance = tol_fixef)
})

test_that("gaussian: fixed-effect priors", {
  prior <- data.frame(
    prior = "normal(0, 0.1)",
    class = "fixef",
    coef = "Days"
  )

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    Reaction ~ Days,
    family = gaussian,
    data = sleepstudy,
    priors = prior,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    Reaction ~ Days,
    family = gaussian,
    data = sleepstudy,
    priors = prior,
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
