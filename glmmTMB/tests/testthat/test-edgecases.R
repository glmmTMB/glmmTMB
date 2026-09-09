stopifnot(require("testthat"),
          require("glmmTMB"))

context("test edge cases")

test_that("profiling failure", {
    ## data from https://github.com/glmmTMB/glmmTMB/issues/399
    dd <- readRDS(system.file("test_data","IC_comp_data.rds", package="glmmTMB"))
    expect_warning(glmmTMB(
        ProbDiv ~ stdQlty + stdLaying + (1|Year) + (1|Site) + (1|PairID),
        family = "binomial",
        control=glmmTMBControl(profile = TRUE),
        data = dd),
        "a Newton step failed")
})

test_that("model with zero total parameters doesn't crash", {
    ## https://github.com/glmmTMB/glmmTMB/issues/1325 : a model with no
    ## fixed effects, no random effects, and a family with no dispersion
    ## parameter has nothing to estimate. TMB/CppAD can't compute a
    ## gradient (or sdreport()) for a zero-length parameter vector, which
    ## used to crash the R process (SIGFPE) rather than fitting trivially.
    d <- data.frame(a = c(5, 3, 7), b = c(5, 7, 3))
    fit <- glmmTMB(cbind(a, b) ~ 0, family = binomial, data = d)
    expect_equal(as.numeric(logLik(fit)), sum(dbinom(d$a, d$a + d$b, 0.5, log = TRUE)))
    expect_length(fixef(fit)$cond, 0)
    expect_equal(unname(predict(fit, type = "response")), rep(0.5, nrow(d)))
    ## downstream methods that rely on vcov()/sdreport() shouldn't crash either
    expect_length(vcov(fit)$cond, 0)
    expect_s3_class(summary(fit), "summary.glmmTMB")

    ## same underlying issue with a family that has no dispersion parameter
    d2 <- data.frame(y = c(5, 3, 7))
    fit2 <- glmmTMB(y ~ 0, family = poisson, data = d2)
    expect_equal(as.numeric(logLik(fit2)), sum(dpois(d2$y, 1, log = TRUE)))

    ## same underlying issue (obj$par has length 0), but reached via
    ## map = fixing every parameter rather than an empty design matrix
    fit3 <- glmmTMB(y ~ 1, family = poisson, data = d2,
                     map = list(beta = factor(NA)), start = list(beta = 0))
    expect_equal(as.numeric(logLik(fit3)), sum(dpois(d2$y, 1, log = TRUE)))
})


