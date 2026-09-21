## Test cases for the RTMB Conway-Maxwell-Binomial family
## Fit each model with RTMB and legacy TMB, then compare likelihoods,
## parameters, predictions, and covariance estimates where applicable.

context("RTMB combinomial backend")

skip_if_not_installed("RTMB")

tol_logLik <- 1e-5
tol_fixef <- 1e-5
tol_varcorr <- 1e-4

rcombinom_test <- function(mean, nu, size) {
  vapply(seq_along(mean), function(i) {
    values <- 0:size[i]
    sample(
      values,
      size = 1L,
      prob = glmmTMB::dcombinom(
        values,
        size = size[i],
        mu = mean[i],
        nu = nu[i]
      )
    )
  }, numeric(1))
}

set.seed(2301)
combinomial_dat <- expand.grid(
  obs = seq_len(5),
  g = factor(seq_len(25))
)
combinomial_dat$x <- rnorm(nrow(combinomial_dat))
combinomial_dat$z <- rnorm(nrow(combinomial_dat))
combinomial_dat$off <- runif(nrow(combinomial_dat), -0.15, 0.15)
combinomial_dat$w <- runif(nrow(combinomial_dat), 0.5, 1.5)
combinomial_dat$size <- sample(8:12, nrow(combinomial_dat), replace = TRUE)
cond_effect <- rnorm(nlevels(combinomial_dat$g), sd = 0.25)
probability <- plogis(
  -0.2 + 0.35 * combinomial_dat$x + cond_effect[combinomial_dat$g]
)
nu <- rep(0.7, nrow(combinomial_dat))
combinomial_dat$success <- rcombinom_test(
  probability * combinomial_dat$size,
  nu,
  combinomial_dat$size
)
combinomial_dat$failure <- combinomial_dat$size - combinomial_dat$success

set.seed(2302)
combinomial_zi_dat <- combinomial_dat
zi_probability <- plogis(-1.8 + 0.3 * combinomial_zi_dat$x)
is_structural_zero <- rbinom(nrow(combinomial_zi_dat), 1L, zi_probability)
combinomial_zi_dat$success[is_structural_zero == 1L] <- 0
combinomial_zi_dat$failure <- combinomial_zi_dat$size -
  combinomial_zi_dat$success

expect_combinomial_matches <- function(m_rtmb, m_tmb,
                                       check_varcorr = FALSE) {
  expect_equal(
    as.numeric(logLik(m_rtmb)),
    as.numeric(logLik(m_tmb)),
    tolerance = tol_logLik
  )
  expect_equal(fixef(m_rtmb)$cond, fixef(m_tmb)$cond, tolerance = tol_fixef)
  expect_equal(fixef(m_rtmb)$disp, fixef(m_tmb)$disp, tolerance = tol_fixef)
  expect_equal(fixef(m_rtmb)$zi, fixef(m_tmb)$zi, tolerance = tol_fixef)
  if (check_varcorr) {
    expect_equal(
      as.numeric(VarCorr(m_rtmb)$cond$g),
      as.numeric(VarCorr(m_tmb)$cond$g),
      tolerance = tol_varcorr
    )
  }
}

test_that("combinomial: RTMB density matches package density", {
  x <- 0:12
  for (nu in c(-0.5, 0.5, 1, 2)) {
    expect_equal(
      glmmTMB:::dcombinom2_rtmb(
        x,
        size = rep(12, length(x)),
        mean = rep(6, length(x)),
        nu = rep(nu, length(x)),
        log = TRUE
      ),
      glmmTMB::dcombinom(x, size = 12, mu = 6, nu = nu, log = TRUE),
      tolerance = 1e-12
    )
  }
})

test_that("combinomial: fixed effects, offsets, and weights", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    cbind(success, failure) ~ x + offset(off),
    weights = w,
    family = combinomial(),
    data = combinomial_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    cbind(success, failure) ~ x + offset(off),
    weights = w,
    family = combinomial(),
    data = combinomial_dat,
    se = FALSE
  )

  expect_combinomial_matches(m_rtmb, m_tmb)
})

test_that("combinomial: dispersion fixed effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    cbind(success, failure) ~ x,
    dispformula = ~ z,
    family = combinomial(),
    data = combinomial_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    cbind(success, failure) ~ x,
    dispformula = ~ z,
    family = combinomial(),
    data = combinomial_dat,
    se = FALSE
  )

  expect_combinomial_matches(m_rtmb, m_tmb)
})

test_that("combinomial: conditional random effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    cbind(success, failure) ~ x + (1 | g),
    family = combinomial(),
    data = combinomial_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    cbind(success, failure) ~ x + (1 | g),
    family = combinomial(),
    data = combinomial_dat,
    se = FALSE
  )

  expect_combinomial_matches(m_rtmb, m_tmb, check_varcorr = TRUE)
})

test_that("combinomial: zero-inflation fixed effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    cbind(success, failure) ~ x,
    ziformula = ~ x,
    family = combinomial(),
    data = combinomial_zi_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    cbind(success, failure) ~ x,
    ziformula = ~ x,
    family = combinomial(),
    data = combinomial_zi_dat,
    se = FALSE
  )

  expect_combinomial_matches(m_rtmb, m_tmb)
})

test_that("combinomial: negative dispersion", {
  start <- list(betadisp = -0.5)
  map <- list(betadisp = factor(NA))

  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    cbind(success, failure) ~ x,
    family = combinomial(allow_negative_nu = TRUE),
    data = combinomial_dat,
    start = start,
    map = map,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    cbind(success, failure) ~ x,
    family = combinomial(allow_negative_nu = TRUE),
    data = combinomial_dat,
    start = start,
    map = map,
    se = FALSE
  )

  expect_combinomial_matches(m_rtmb, m_tmb)
  expect_equal(sigma(m_rtmb), -0.5, tolerance = tol_fixef)
})

test_that("combinomial: prediction and simulation", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    cbind(success, failure) ~ x,
    family = combinomial(),
    data = combinomial_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    cbind(success, failure) ~ x,
    family = combinomial(),
    data = combinomial_dat,
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

  set.seed(2303)
  simulated <- simulate(m_rtmb, nsim = 1)[[1]]
  expect_equal(dim(simulated), c(nrow(combinomial_dat), 2L))
  expect_equal(rowSums(simulated), combinomial_dat$size)
})
