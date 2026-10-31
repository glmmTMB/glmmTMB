## Tests for RTMB random-effect simulation controls

context("RTMB random-effect simulation controls")

skip_if_not_installed("RTMB")

data("sleepstudy", package = "lme4")

test_that("random-effect simCode modes match TMB", {
  for (backend in c(TRUE, FALSE)) {
    glmmTMB::useRTMB(backend)
    model <- glmmTMB(
      Reaction ~ Days + (1 | Subject),
      family = gaussian,
      data = sleepstudy,
      se = FALSE
    )
    fitted_b <- model$obj$report()$b

    set_simcodes(model$obj, "zero")
    zero_simulation <- model$obj$simulate(complete = TRUE)
    expect_equal(
      zero_simulation$b,
      rep(0, length(fitted_b)),
      info = paste("backend:", backend)
    )

    set_simcodes(model$obj, "fix")
    fixed_simulation <- model$obj$simulate(complete = TRUE)
    expect_equal(
      fixed_simulation$b,
      fitted_b,
      info = paste("backend:", backend)
    )

    set_simcodes(model$obj, "random")
    random_simulation <- model$obj$simulate(complete = TRUE)
    expect_true(
      all(is.finite(random_simulation$b)),
      info = paste("backend:", backend)
    )
    expect_false(
      isTRUE(all.equal(random_simulation$b, fitted_b)),
      info = paste("backend:", backend)
    )
  }
})

test_that("random-only covariance structures reject other simCode modes", {
  simulation_data <- transform(sleepstudy, DaysFac = factor(Days))

  for (backend in c(TRUE, FALSE)) {
    glmmTMB::useRTMB(backend)
    model <- glmmTMB(
      Reaction ~ Days + homdiag(0 + DaysFac | Subject),
      family = gaussian,
      data = simulation_data,
      se = FALSE
    )

    set_simcodes(model$obj, "zero")
    expect_error(
      model$obj$simulate(complete = TRUE),
      "sim[Cc]ode.*implemented",
      info = paste("backend:", backend)
    )
  }
})

test_that("simulate() works with REML = TRUE and matches TMB (GH #1358)", {
  sims <- lapply(c(rtmb = TRUE, tmb = FALSE), function(backend) {
    glmmTMB::useRTMB(backend)
    model <- glmmTMB(
      Reaction ~ Days + (1 | Subject),
      data = sleepstudy,
      REML = TRUE,
      se = FALSE
    )
    simulate(model, nsim = 2, seed = 5006)
  })
  expect_equal(sims$rtmb, sims$tmb)
})
