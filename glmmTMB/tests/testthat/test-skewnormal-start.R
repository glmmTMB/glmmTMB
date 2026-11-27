test_that("zero-start skewnormal fits can find either sign (GH #1375)", {
    for (shape in c(-3, 3)) {
        dd <- data.frame(id = seq_len(2000))
        dd$y <- simulate_new(~1, family = skewnormal, newdata = dd,
                            newparams = list(beta = 1, betadisp = log(2),
                                             psi = shape), seed = 1)[[1]]
        fit <- glmmTMB(y ~ 1, family = skewnormal, data = dd)
        reference <- update(fit, start = list(psi = shape))
        symmetric <- update(fit, map = list(psi = factor(NA)))

        expect_gt(shape * family_params(fit), 1)
        expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(reference)),
                     tolerance = 1e-6)
        expect_gt(as.numeric(logLik(fit) - logLik(symmetric)), 10)
        expect_true(fit$sdr$pdHess)
        expect_equal(fit$fit$parfull, fit$fit$par)
        expect_equal(unname(residuals(fit, type = "pearson")),
                     unname((dd$y - fitted(fit)) / sigma(fit)))
    }
    optimizer <- function(...) {
        res <- nlminb(...)
        res$par <- unname(res$par)
        res
    }
    unnamed <- update(fit, control = glmmTMBControl(optimizer = optimizer))
    expect_equal(as.numeric(logLik(unnamed)), as.numeric(logLik(fit)),
                 tolerance = 1e-6)
})

test_that("skewnormal initialization follows residual rather than response skew", {
    set.seed(123)
    dd <- data.frame(x = rgamma(2000, 2))
    dd$y <- simulate_new(~x, family = skewnormal, newdata = dd,
                        newparams = list(beta = c(1, 2), betadisp = log(2),
                                         psi = -3), seed = 1)[[1]]
    expect_gt(mean((dd$y - mean(dd$y))^3), 0)
    fit <- glmmTMB(y ~ x, family = skewnormal, data = dd)
    reference <- update(fit, start = list(psi = -3))
    expect_lt(family_params(fit), -1)
    expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(reference)),
                 tolerance = 1e-6)
    expect_true(fit$sdr$pdHess)
})

test_that("skewnormal restarts respect supplied starts and fixed shape", {
    dd <- data.frame(y = seq(-2, 2, length.out = 100))
    starts <- list()
    optimizer <- function(par, ...) {
        starts[[length(starts) + 1L]] <<- par
        nlminb(par, ...)
    }
    ctrl <- glmmTMBControl(optimizer = optimizer)
    glmmTMB(y ~ 1, family = skewnormal, data = dd, control = ctrl)
    expect_equal(vapply(starts, function(p) unname(p["psi"]), numeric(1)),
                 c(0, -1, 1))

    starts <- list()
    glmmTMB(y ~ 1, family = skewnormal, data = dd, control = ctrl,
            start = list(psi = 1))
    expect_length(starts, 1L)
    expect_equal(unname(starts[[1]]["psi"]), 1)

    starts <- list()
    fit <- glmmTMB(y ~ 1, family = skewnormal, data = dd, control = ctrl,
                   map = list(psi = factor(NA)))
    expect_length(starts, 1L)
    expect_equal(unname(family_params(fit)), 0)

    starts <- list()
    glmmTMB(y ~ 1, data = dd, control = ctrl)
    expect_length(starts, 1L)
})

test_that("skewnormal restarts work with profiling and RTMB", {
    dd <- data.frame(id = seq_len(1000))
    dd$y <- simulate_new(~1, family = skewnormal, newdata = dd,
                        newparams = list(beta = 1, betadisp = 0,
                                         psi = 3), seed = 1)[[1]]
    fit <- glmmTMB(y ~ 1, family = skewnormal, data = dd)
    profiled <- update(fit, control = glmmTMBControl(profile = TRUE))
    expect_equal(as.numeric(logLik(profiled)), as.numeric(logLik(fit)),
                 tolerance = 1e-6)
    expect_equal(profiled$fit$parfull, profiled$fit$par)

    local_useRTMB(TRUE)
    rtmb <- update(fit)
    expect_equal(as.numeric(logLik(rtmb)), as.numeric(logLik(fit)),
                 tolerance = 1e-6)
    expect_equal(family_params(rtmb), family_params(fit), tolerance = 1e-4)
})

test_that("a failed skewnormal restart does not replace the retained fit", {
    dd <- data.frame(id = seq_len(1000))
    dd$y <- simulate_new(~1, family = skewnormal, newdata = dd,
                        newparams = list(beta = 1, betadisp = 0,
                                         psi = 3), seed = 1)[[1]]
    optimizer <- function(par, ...) {
        res <- nlminb(par, ...)
        if (par["psi"] > 0) stop("failed restart")
        res
    }
    fit <- glmmTMB(y ~ 1, family = skewnormal, data = dd, se = FALSE,
                   control = glmmTMBControl(optimizer = optimizer))
    expect_equal(fit$fit$parfull, fit$fit$par)
    expect_equal(fit$obj$fn(fit$fit$par), fit$fit$objective)
})

test_that("skewnormal restarts retain the selected random effects", {
    set.seed(123)
    dd <- expand.grid(id = seq_len(50), g = factor(seq_len(20)))
    dd$x <- rnorm(nrow(dd))
    dd$y <- simulate_new(~x + (1 | g), family = skewnormal, newdata = dd,
                        newparams = list(beta = c(1, 0.3), betadisp = 0,
                                         theta = log(0.5), psi = 3),
                        seed = 1)[[1]]
    fit <- glmmTMB(y ~ x + (1 | g), family = skewnormal, data = dd)
    reference <- update(fit, start = list(psi = 3))
    expect_gt(family_params(fit), 1)
    expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(reference)),
                 tolerance = 1e-6)
    expect_equal(fit$fit$parfull[-fit$obj$env$random], fit$fit$par)
    expect_equal(getME(fit, "b"), getME(reference, "b"), tolerance = 1e-4)
    expect_true(fit$sdr$pdHess)
})
