stopifnot(require("testthat"),
          require("glmmTMB"))

data(sleepstudy, cbpp, Pastes,
     package = "lme4")

test_that("diagnose works with REML fits", {
    fm1 <- glmmTMB(Reaction ~ Days + (1|Subject),
                   sleepstudy,
                   REML = TRUE)
    cc <- capture.output(d <- diagnose(fm1, explain = FALSE))
    expect_false(d)
})

test_that("diagnose works with Tweedie fits", {
    skip_on_cran()
    ## false convergence
    fm1 <- suppressWarnings(glmmTMB(Reaction ~ Days + (1|Subject),
                   sleepstudy,
                   family = tweedie))
    cc <- capture.output(d <- diagnose(fm1, explain = FALSE))
    expect_false(d)
})
