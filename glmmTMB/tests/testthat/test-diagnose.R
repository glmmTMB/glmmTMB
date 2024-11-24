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

test_that("diagnose with mapped parameters", {
    ## GH 1120
    M0 <- suppressWarnings(
        glmmTMB(
            Reaction ~ Days + cs(0 + factor(Days) | Subject), dispformula = ~ 0,
            map = list(theta = factor(c(rep(1, 10), 2))),
            data = sleepstudy, REML = TRUE)
    )
    dd_out <- capture.output(dd <- diagnose(M0))
    expect_equal(dd, FALSE)
    expect_true(any(grepl("Unusually large Z-statistics", dd_out)))
})
