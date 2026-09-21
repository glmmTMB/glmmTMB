## Test cases for the RTMB betabinomial family
## Fit each model with the RTMB and legacy TMB backends, then compare
## likelihoods, fixed effects, predictions, and covariance estimates where
## applicable.

context("RTMB betabinomial backend")

skip_if_not_installed("RTMB")

tol_logLik <- 1e-5
tol_fixef <- 1e-5
tol_varcorr <- 1e-4

set.seed(1901)
bb_dat <- expand.grid(
  obs = seq_len(5),
  g = factor(seq_len(30))
)
bb_dat$x <- rnorm(nrow(bb_dat))
bb_dat$off <- runif(nrow(bb_dat), -0.2, 0.2)
bb_dat$size <- 10
bb_effect <- rnorm(nlevels(bb_dat$g), sd = 0.3)
bb_prob <- plogis(-0.3 + 0.35 * bb_dat$x + bb_effect[bb_dat$g])
bb_phi <- 7
bb_draw_prob <- rbeta(
  nrow(bb_dat),
  shape1 = bb_prob * bb_phi,
  shape2 = (1 - bb_prob) * bb_phi
)
bb_dat$y <- rbinom(nrow(bb_dat), size = bb_dat$size, prob = bb_draw_prob)

set.seed(1902)
bb_zi_dat <- bb_dat
zi_effect <- rnorm(nlevels(bb_zi_dat$g), sd = 0.4)
zi_prob <- plogis(-1.4 + zi_effect[bb_zi_dat$g])
is_structural_zero <- rbinom(nrow(bb_zi_dat), size = 1, prob = zi_prob)
bb_zi_dat$y <- ifelse(is_structural_zero == 1, 0, bb_zi_dat$y)

expect_betabinomial_equal <- function(m_rtmb, m_tmb, check_varcorr = FALSE) {
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
    fixef(m_rtmb)$disp,
    fixef(m_tmb)$disp,
    tolerance = tol_fixef
  )
  expect_equal(
    fixef(m_rtmb)$zi,
    fixef(m_tmb)$zi,
    tolerance = tol_fixef
  )
  if (check_varcorr) {
    expect_equal(
      unname(unlist(VarCorr(m_rtmb)$cond)),
      unname(unlist(VarCorr(m_tmb)$cond)),
      tolerance = tol_varcorr
    )
  }
}

test_that("betabinomial: RTMB density matches beta-binomial formula", {
  x <- 0:5
  size <- rep(5, length(x))
  shape1 <- 2.5
  shape2 <- 4

  expect_equal(
    glmmTMB:::dbetabinom_robust_rtmb(
      x,
      log_shape1 = log(shape1),
      log_shape2 = log(shape2),
      size = size,
      log = TRUE
    ),
    lchoose(size, x) + lbeta(x + shape1, size - x + shape2) -
      lbeta(shape1, shape2),
    tolerance = 1e-12
  )
})

test_that("betabinomial: grouped response with weights", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y / size ~ x,
    weights = size,
    family = betabinomial(link = "logit"),
    data = bb_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y / size ~ x,
    weights = size,
    family = betabinomial(link = "logit"),
    data = bb_dat,
    se = FALSE
  )

  expect_betabinomial_equal(m_rtmb, m_tmb)
})

test_that("betabinomial: cbind response", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    cbind(y, size - y) ~ x,
    family = betabinomial(link = "logit"),
    data = bb_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    cbind(y, size - y) ~ x,
    family = betabinomial(link = "logit"),
    data = bb_dat,
    se = FALSE
  )

  expect_betabinomial_equal(m_rtmb, m_tmb)
})

test_that("betabinomial: offsets", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y / size ~ x + offset(off),
    weights = size,
    family = betabinomial(link = "logit"),
    data = bb_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y / size ~ x + offset(off),
    weights = size,
    family = betabinomial(link = "logit"),
    data = bb_dat,
    se = FALSE
  )

  expect_betabinomial_equal(m_rtmb, m_tmb)
})

test_that("betabinomial: dispersion fixed effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y / size ~ x,
    weights = size,
    dispformula = ~ x,
    family = betabinomial(link = "logit"),
    data = bb_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y / size ~ x,
    weights = size,
    dispformula = ~ x,
    family = betabinomial(link = "logit"),
    data = bb_dat,
    se = FALSE
  )

  expect_betabinomial_equal(m_rtmb, m_tmb)
})

test_that("betabinomial: conditional random intercept", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y / size ~ x + (1 | g),
    weights = size,
    family = betabinomial(link = "logit"),
    data = bb_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y / size ~ x + (1 | g),
    weights = size,
    family = betabinomial(link = "logit"),
    data = bb_dat,
    se = FALSE
  )

  expect_betabinomial_equal(m_rtmb, m_tmb, check_varcorr = TRUE)
})

test_that("betabinomial: cloglog link", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y / size ~ x,
    weights = size,
    family = betabinomial(link = "cloglog"),
    data = bb_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y / size ~ x,
    weights = size,
    family = betabinomial(link = "cloglog"),
    data = bb_dat,
    se = FALSE
  )

  expect_betabinomial_equal(m_rtmb, m_tmb)
})

test_that("betabinomial: zero-inflation fixed effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y / size ~ x,
    weights = size,
    ziformula = ~ x,
    family = betabinomial(link = "logit"),
    data = bb_zi_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y / size ~ x,
    weights = size,
    ziformula = ~ x,
    family = betabinomial(link = "logit"),
    data = bb_zi_dat,
    se = FALSE
  )

  expect_betabinomial_equal(m_rtmb, m_tmb)
})

test_that("betabinomial: zero-inflation random effects", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y / size ~ x,
    weights = size,
    ziformula = ~ 0 + (1 | g),
    family = betabinomial(link = "logit"),
    data = bb_zi_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y / size ~ x,
    weights = size,
    ziformula = ~ 0 + (1 | g),
    family = betabinomial(link = "logit"),
    data = bb_zi_dat,
    se = FALSE
  )

  expect_betabinomial_equal(m_rtmb, m_tmb)
})

test_that("betabinomial: prediction and simulation", {
  local_useRTMB(TRUE)
  m_rtmb <- glmmTMB(
    y / size ~ x,
    weights = size,
    family = betabinomial(link = "logit"),
    data = bb_dat,
    se = FALSE
  )

  glmmTMB::useRTMB(FALSE)
  m_tmb <- glmmTMB(
    y / size ~ x,
    weights = size,
    family = betabinomial(link = "logit"),
    data = bb_dat,
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

  set.seed(1903)
  sim_y <- simulate(m_rtmb, nsim = 1)[[1]]
  expect_equal(nrow(sim_y), nrow(bb_dat))
  expect_equal(rowSums(sim_y), bb_dat$size)
  expect_true(all(sim_y[, 1] >= 0 & sim_y[, 1] <= bb_dat$size))
  expect_true(all(sim_y[, 1] == floor(sim_y[, 1])))
})
