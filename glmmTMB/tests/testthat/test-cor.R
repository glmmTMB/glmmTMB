## Tests for the 'cor' covariance structure (unstructured correlation,
## unit variances; GH #1356)

stopifnot(require("testthat"),
          require("glmmTMB"))

set.seed(20261031)
ng <- 400
nt <- 3
dd_cor <- expand.grid(time = factor(seq_len(nt)), id = factor(seq_len(ng)))
dd_cor$y <- 0
Cmat <- matrix(c(1, 0.6, 0.3,
                 0.6, 1, 0.5,
                 0.3, 0.5, 1), nt)
pars_cor <- list(beta = c(0, 0.5, -0.5), theta = put_cor(Cmat))
form_cor <- ~ time + cor(0 + time | id)
dd_cor$y <- simulate_new(form_cor, newdata = dd_cor,
                         family = binomial, weights = rep(20, nrow(dd_cor)),
                         newparams = pars_cor, seed = 101)[[1]] / 20

fit_cor <- glmmTMB(y ~ time + cor(0 + time | id), data = dd_cor,
                   family = binomial, weights = rep(20, nrow(dd_cor)))

test_that("cor structure has n(n-1)/2 parameters and unit SDs", {
    expect_equal(length(getME(fit_cor, "theta")), nt * (nt - 1) / 2)
    vc <- VarCorr(fit_cor)$cond$id
    expect_equal(unname(attr(vc, "stddev")), rep(1, nt))
    expect_equal(unname(diag(vc)), rep(1, nt))
})

test_that("cor matches us with log-SDs mapped to zero", {
    fit_us <- glmmTMB(y ~ time + us(0 + time | id), data = dd_cor,
                      family = binomial, weights = rep(20, nrow(dd_cor)),
                      start = list(theta = rep(0, 6)),
                      map = list(theta = factor(c(NA, NA, NA, 1:3))))
    expect_equal(c(logLik(fit_cor)), c(logLik(fit_us)), tolerance = 1e-6)
    expect_equal(fixef(fit_cor)$cond, fixef(fit_us)$cond, tolerance = 1e-5)
    expect_equal(unname(getME(fit_cor, "theta")),
                 unname(getME(fit_us, "theta"))[4:6], tolerance = 1e-4)
    vc_cor <- VarCorr(fit_cor)$cond$id
    vc_us <- VarCorr(fit_us)$cond$id
    expect_equal(c(vc_cor), c(vc_us), tolerance = 1e-4)
    expect_equal(attr(vc_cor, "correlation"), attr(vc_us, "correlation"),
                 tolerance = 1e-4)
})

test_that("cor fit recovers simulated correlations", {
    cc <- attr(VarCorr(fit_cor)$cond$id, "correlation")
    expect_equal(cc[lower.tri(cc)], Cmat[lower.tri(Cmat)], tolerance = 0.1)
})

test_that("cor structure requires blocks of size > 1", {
    expect_error(glmmTMB(y ~ time + cor(1 | id), data = dd_cor,
                         family = binomial, weights = rep(20, nrow(dd_cor))),
                 "requires at least two")
})

test_that("simulate_new with cor gives unit-variance correlated REs", {
    obj <- simulate_new(form_cor, newdata = dd_cor, family = binomial,
                        newparams = pars_cor, return_val = "object")
    vc <- VarCorr(obj)$cond$id
    expect_equal(unname(attr(vc, "stddev")), rep(1, nt))
    expect_equal(unname(attr(vc, "correlation")), Cmat, tolerance = 1e-8)
    ## Gaussian response with negligible residual variance: response = b
    b <- simulate_new(form_cor, newdata = dd_cor, family = gaussian,
                      newparams = c(list(betadisp = -10),
                                    replace(pars_cor, "beta", list(rep(0, nt)))),
                      seed = 102)[[1]]
    bmat <- matrix(b, ncol = nt, byrow = TRUE)
    expect_equal(cov(bmat), Cmat, tolerance = 0.15)
})

test_that("cor structure supports all simCode modes", {
    fitted_b <- fit_cor$obj$report()$b
    set_simcodes(fit_cor$obj, "zero")
    expect_equal(fit_cor$obj$simulate(complete = TRUE)$b,
                 rep(0, length(fitted_b)))
    set_simcodes(fit_cor$obj, "fix")
    expect_equal(fit_cor$obj$simulate(complete = TRUE)$b, fitted_b)
    set_simcodes(fit_cor$obj, "random")
    bsim <- fit_cor$obj$simulate(complete = TRUE)$b
    expect_true(all(is.finite(bsim)))
    expect_false(isTRUE(all.equal(bsim, fitted_b)))
    ss <- simulate(fit_cor, nsim = 2, seed = 103)
    expect_equal(dim(ss), c(nrow(dd_cor), 2L))
})

test_that("methods work for cor fits", {
    expect_output(print(VarCorr(fit_cor)), "Corr")
    expect_s3_class(summary(fit_cor), "summary.glmmTMB")
    ci <- confint(fit_cor)
    ## SDs are fixed, so only correlation CIs remain
    expect_equal(sum(grepl("Cor\\.", rownames(ci))), 3L)
    expect_false(any(grepl("Std\\.Dev", rownames(ci))))
    expect_length(predict(fit_cor), nrow(dd_cor))
})

test_that("ranef_sd priors are rejected for cor terms", {
    pr <- data.frame(prior = "normal(0, 1)", class = "ranef_sd", coef = "id")
    expect_error(update(fit_cor, priors = pr), "no SD parameters")
    pr2 <- data.frame(prior = "lkj(2)", class = "ranef_cor", coef = "id")
    expect_no_error(update(fit_cor, priors = pr2))
    ## blanket SD priors (blank or omitted coef) must not reach correlations
    pr_blank <- data.frame(prior = "normal(0, 1)", class = "ranef_sd", coef = "")
    expect_error(update(fit_cor, priors = pr_blank),
                 "no random-effects terms.*standard-deviation")
    pr_nocoef <- data.frame(prior = "normal(0, 1)", class = "ranef_sd")
    expect_error(update(fit_cor, priors = pr_nocoef),
                 "no random-effects terms.*standard-deviation")
    ## a blanket correlation prior covers all three correlation parameters
    pr_cor <- data.frame(prior = "lkj(2)", class = "ranef_cor")
    fit_pr <- update(fit_cor, priors = pr_cor)
    dd <- fit_pr$obj$env$data
    expect_equal(c(dd$prior_elstart, dd$prior_elend), c(0, 2))
})

test_that("cor respects full_cor = FALSE on both back ends", {
    for (backend in c(FALSE, TRUE)) {
        if (backend) skip_if_not_installed("RTMB")
        old <- glmmTMB::useRTMB()
        glmmTMB::useRTMB(backend)
        fit <- glmmTMB(y ~ time + cor(0 + time | id), data = dd_cor,
                       family = binomial, weights = rep(20, nrow(dd_cor)),
                       se = FALSE, control = glmmTMBControl(full_cor = FALSE))
        glmmTMB::useRTMB(old)
        expect_identical(fit$obj$report()$corr[[1]], matrix(NaN, 1, 1),
                         info = paste("RTMB:", backend))
    }
})

test_that("RTMB backend matches TMB for cor", {
    skip_if_not_installed("RTMB")
    old <- glmmTMB::useRTMB()
    on.exit(glmmTMB::useRTMB(old))
    glmmTMB::useRTMB(TRUE)
    fit_rtmb <- glmmTMB(y ~ time + cor(0 + time | id), data = dd_cor,
                        family = binomial, weights = rep(20, nrow(dd_cor)),
                        se = FALSE)
    expect_equal(c(logLik(fit_rtmb)), c(logLik(fit_cor)), tolerance = 1e-5)
    expect_equal(fixef(fit_rtmb)$cond, fixef(fit_cor)$cond, tolerance = 1e-5)
    expect_equal(VarCorr(fit_rtmb), VarCorr(fit_cor), tolerance = 1e-4)
    fitted_b <- fit_rtmb$obj$report()$b
    set_simcodes(fit_rtmb$obj, "zero")
    expect_equal(fit_rtmb$obj$simulate(complete = TRUE)$b,
                 rep(0, length(fitted_b)))
    set_simcodes(fit_rtmb$obj, "random")
    bsim <- fit_rtmb$obj$simulate(complete = TRUE)$b
    expect_equal(length(bsim), length(fitted_b))
    expect_true(all(is.finite(bsim)))
})
