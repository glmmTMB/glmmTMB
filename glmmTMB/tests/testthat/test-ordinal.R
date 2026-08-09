## ordinal (cumulative link / proportional odds) family

stopifnot(require("testthat"),
          require("glmmTMB"))

data("housing", package = "MASS")

fit_ord <- glmmTMB(Sat ~ Infl + Type + Cont, weights = Freq,
                   data = housing, family = ordinal())

test_that("ordinal fixed effects match MASS::polr", {
    fit_polr <- MASS::polr(Sat ~ Infl + Type + Cont, weights = Freq,
                           data = housing)
    expect_equal(unname(fixef(fit_ord)$cond[-1]),
                 unname(coef(fit_polr)),
                 tolerance = 1e-4)
    expect_equal(unname(family_params(fit_ord)),
                 unname(fit_polr$zeta),
                 tolerance = 1e-4)
    expect_equal(c(logLik(fit_ord)), c(logLik(fit_polr)),
                 tolerance = 1e-6)
    ## threshold names built from response levels
    expect_identical(names(family_params(fit_ord)),
                     c("Low|Medium", "Medium|High"))
    ## intercept fixed to zero (absorbed into thresholds)
    expect_equal(unname(fixef(fit_ord)$cond["(Intercept)"]), 0)
    ## ... and suppressed from the summary coefficient table
    expect_false("(Intercept)" %in%
                 rownames(summary(fit_ord)$coefficients$cond))
    ## internal map does not leak into user-facing modelInfo$map
    expect_null(fit_ord$modelInfo$map)
})

test_that("ordinal probit link matches MASS::polr", {
    fit_polr_p <- MASS::polr(Sat ~ Infl + Type + Cont, weights = Freq,
                             data = housing, method = "probit")
    fit_ord_p <- update(fit_ord, family = ordinal(link = "probit"))
    expect_equal(unname(fixef(fit_ord_p)$cond[-1]),
                 unname(coef(fit_polr_p)),
                 tolerance = 1e-4)
    expect_equal(c(logLik(fit_ord_p)), c(logLik(fit_polr_p)),
                 tolerance = 1e-6)
})

test_that("ordinal predictions: probs, response, link", {
    pr <- predict(fit_ord, type = "probs")
    expect_true(is.matrix(pr))
    expect_identical(colnames(pr), levels(housing$Sat))
    expect_equal(unname(rowSums(pr)), rep(1, nrow(pr)))
    ## response prediction is the expected category index
    expect_equal(predict(fit_ord, type = "response"),
                 unname(pr %*% seq_len(ncol(pr)))[, 1])
    ## link prediction is the latent-scale linear predictor:
    ## P(Y <= j) = linkinv(theta_j - eta)
    eta <- predict(fit_ord, type = "link")
    theta <- family_params(fit_ord)
    expect_equal(unname(pr[, 1]), unname(plogis(theta[1] - eta)))
    ## se.fit works and returns matching matrices
    prs <- predict(fit_ord, type = "probs", se.fit = TRUE)
    expect_identical(dim(prs$fit), dim(prs$se.fit))
    expect_false(anyNA(prs$se.fit))
    ## type = "probs" is ordinal-only
    fit_pois <- glmmTMB(count ~ mined, family = poisson, data = Salamanders)
    expect_error(predict(fit_pois, type = "probs"), "only available")
})

test_that("ordinal mixed model matches ordinal::clmm", {
    skip_if_not_installed("ordinal")
    data("wine", package = "ordinal")
    fit_clmm <- ordinal::clmm(rating ~ temp + contact + (1 | judge),
                              data = wine)
    fit_tmb <- glmmTMB(rating ~ temp + contact + (1 | judge),
                       data = wine, family = ordinal())
    expect_equal(c(logLik(fit_tmb)), c(logLik(fit_clmm)), tolerance = 1e-5)
    expect_equal(unname(fixef(fit_tmb)$cond[-1]),
                 unname(fit_clmm$beta), tolerance = 1e-3)
    expect_equal(unname(family_params(fit_tmb)),
                 unname(fit_clmm$alpha), tolerance = 1e-3)
    expect_equal(unname(attr(VarCorr(fit_tmb)$cond$judge, "stddev")),
                 unname(sqrt(ordinal::VarCorr(fit_clmm)$judge[1, 1])),
                 tolerance = 1e-3)
})

test_that("ordinal simulate/residuals/refit", {
    set.seed(101)
    dd <- data.frame(x = rnorm(600),
                     g = factor(rep(1:30, each = 20)))
    eta <- dd$x + rnorm(30)[as.integer(dd$g)]
    u <- runif(600)
    cum <- plogis(outer(c(-1, 0.5, 2), eta, "-"))
    dd$y <- ordered(1 + colSums(sweep(cum, 2, u, "<")), levels = 1:4)
    fit <- glmmTMB(y ~ x + (1 | g), data = dd, family = ordinal())

    ## simulate returns ordered factors on the original levels
    ss <- simulate(fit, nsim = 2, seed = 1)
    expect_true(is.ordered(ss[[1]]))
    expect_identical(levels(ss[[1]]), levels(dd$y))

    ## refit to simulated data works, integer codes accepted as response
    dd2 <- dd
    dd2$y <- as.numeric(ss[[1]])
    fit2 <- update(fit, data = dd2)
    expect_true(fit2$fit$convergence == 0)

    ## Dunn-Smyth residuals approximately standard normal
    set.seed(1)
    r <- residuals(fit, type = "dunn-smyth")
    expect_true(abs(mean(r)) < 0.15 && abs(sd(r) - 1) < 0.15)

    ## response residuals on category-index scale
    rr <- residuals(fit, type = "response")
    expect_equal(unname(rr),
                 as.numeric(dd$y) - predict(fit, type = "response"))
})

test_that("ordinal downstream: emmeans and car::Anova handle mapped intercept", {
    skip_if_not_installed("emmeans")
    em <- emmeans::emmeans(fit_ord, ~ Infl)
    es <- summary(em)
    expect_false(anyNA(es$SE))
    ## contrasts reproduce the fixed-effect coefficients
    ec <- summary(emmeans::contrast(em, "trt.vs.ctrl"))
    expect_equal(ec$estimate,
                 unname(fixef(fit_ord)$cond[c("InflMedium", "InflHigh")]),
                 tolerance = 1e-6)

    skip_if_not_installed("car")
    aa <- car::Anova(fit_ord)
    expect_false(anyNA(aa[["Chisq"]]))
    ## Wald chisq consistent with squared z for the 1-df term
    z_cont <- summary(fit_ord)$coefficients$cond["ContHigh", "z value"]
    expect_equal(aa["Cont", "Chisq"], z_cont^2, tolerance = 1e-6)
})

test_that("ordinal numerical robustness at extreme parameters", {
    ## probit/cloglog cumulative log-probs must stay finite at extreme
    ## eta (logspace_sub(-Inf, -Inf) = NaN would poison the gradient)
    set.seed(1)
    dd <- data.frame(x = c(rnorm(299), 12))
    u <- runif(300)
    cum <- plogis(outer(qlogis(c(.25, .5, .75)), 1.5 * dd$x, "-"))
    dd$y <- ordered(1 + colSums(sweep(cum, 2, u, "<")), levels = 1:4)
    for (lnk in c("probit", "cloglog")) {
        obj <- glmmTMB(y ~ x, data = dd, family = ordinal(link = lnk),
                       doFit = FALSE)
        ff <- fitTMB(obj, doOptim = FALSE)
        p0 <- ff$par
        p0[names(p0) == "beta"] <- 8
        expect_true(is.finite(ff$fn(p0)))
        expect_false(anyNA(ff$gr(p0)))
    }
    ## threshold transform exact under dominating category weights
    f4 <- glmmTMB(y ~ x, data = dd, family = ordinal(),
                  start = list(psi = c(45, 0, 0)), doFit = FALSE)
    ff4 <- fitTMB(f4, doOptim = FALSE)
    expect_true(is.finite(ff4$fn(ff4$par)))
})

test_that("ordinal Anova type III and confint thresholds", {
    skip_if_not_installed("car")
    a3 <- car::Anova(fit_ord, type = 3)
    ## fixed-to-zero intercept is untestable -> NA row, others finite
    expect_true(is.na(a3["(Intercept)", "Chisq"]))
    expect_false(anyNA(a3[c("Infl", "Type", "Cont"), "Chisq"]))
    ## confint includes delta-method threshold CIs
    ci <- confint(fit_ord, component = "all")
    expect_true(all(c("Low|Medium", "Medium|High") %in% rownames(ci)))
    thr <- family_params(fit_ord)
    expect_true(all(ci[names(thr), 1] < thr & thr < ci[names(thr), 2]))
})

test_that("ordinal integer-coded responses: warning and numeric simulate", {
    dd <- data.frame(x = rnorm(300))
    set.seed(3)
    dd$y <- 1 + rbinom(300, 3, plogis(dd$x))
    expect_warning(fit_i <- glmmTMB(y ~ x, data = dd, family = ordinal()),
                   "integer codes")
    ## simulate preserves the numeric response type
    si <- simulate(fit_i, nsim = 1, seed = 1)[[1]]
    expect_true(is.numeric(si))
    expect_true(all(si %in% 1:4))
})

test_that("ordinal error handling", {
    expect_error(glmmTMB(Sat ~ Infl, weights = Freq, data = housing,
                         ziformula = ~1, family = ordinal()),
                 "zero-inflation is not implemented")
    housing2 <- transform(housing, Sat0 = as.numeric(Sat) - 1)
    expect_error(glmmTMB(Sat0 ~ Infl, weights = Freq, data = housing2,
                         family = ordinal()),
                 "ordered factor")
    housing3 <- transform(housing, Satu = factor(as.character(Sat)))
    expect_warning(glmmTMB(Satu ~ Infl, weights = Freq, data = housing3,
                           family = ordinal()),
                   "unordered factor")
})

test_that("ordinal REML warns that thresholds are not integrated out", {
    expect_warning(fit <- glmmTMB(Sat ~ Infl, weights = Freq, data = housing,
                                  family = ordinal(), REML = TRUE),
                   "not the thresholds")
    ## the fit still succeeds: this is a caveat, not a prohibition
    expect_s3_class(fit, "glmmTMB")
    expect_true(fit$modelInfo$REML)
    ## ML is unaffected
    expect_silent(glmmTMB(Sat ~ Infl, weights = Freq, data = housing,
                          family = ordinal()))
})
