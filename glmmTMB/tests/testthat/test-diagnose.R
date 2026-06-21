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

test_that("diagnose scale report shows the threshold, not the flagged SDs", {
    ## the check_scales header should interpolate big_sd_log10 (the threshold),
    ## not the vector of flagged standard deviations
    set.seed(101)
    n <- 300
    dd <- data.frame(y       = rnorm(n),
                     x_big   = rnorm(n, sd = 1e4),    # |log10(sd)| ~ 4 > 3
                     x_small = rnorm(n, sd = 1e-4))   # |log10(sd)| ~ 4 > 3
    fm <- suppressWarnings(glmmTMB(y ~ x_big + x_small, data = dd))
    out <- suppressWarnings(capture.output(
        diagnose(fm, explain = FALSE,
                 check_coefs = FALSE, check_zstats = FALSE,
                 check_hessian = FALSE, check_scales = TRUE)))
    hdr <- grep("standard deviations", out, value = TRUE)
    ## exactly one header line (the bug recycled it once per flagged predictor) ...
    expect_length(hdr, 1L)
    ## ... carrying the threshold big_sd_log10 = 3, not a flagged SD value
    expect_match(hdr, "|log10(sd)|>3):", fixed = TRUE)
})
