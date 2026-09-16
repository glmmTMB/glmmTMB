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
    ## threshold SEs on the theta scale match clmm (issue #1323)
    expect_equal(unname(glmmTMB:::ordinal_thresholds(fit_tmb)[, "Std. Error"]),
                 unname(sqrt(diag(vcov(fit_clmm)))[seq_along(fit_clmm$alpha)]),
                 tolerance = 1e-3)
})

test_that("ordinal threshold standard errors (delta method)", {
    fit_polr <- MASS::polr(Sat ~ Infl + Type + Cont, weights = Freq,
                           data = housing, Hess = TRUE)
    thr <- glmmTMB:::ordinal_thresholds(fit_ord)
    expect_identical(rownames(thr), c("Low|Medium", "Medium|High"))
    expect_equal(thr[, "Estimate"], family_params(fit_ord))
    expect_equal(unname(thr[, "Std. Error"]),
                 unname(summary(fit_polr)$coefficients[rownames(thr),
                                                       "Std. Error"]),
                 tolerance = 1e-4)
    ## consistent with the Wald CIs from confint()
    ci <- confint(fit_ord, component = "all")
    expect_equal(unname(thr[, "Std. Error"]),
                 unname((ci[rownames(thr), 2] - ci[rownames(thr), 1]) /
                        (2 * qnorm(0.975))),
                 tolerance = 1e-8)
    ## exposed in summary() as a separate 'thresholds' table
    ss <- summary(fit_ord)
    expect_identical(colnames(ss$thresholds),
                     c("Estimate", "Std. Error", "z value"))
    expect_equal(ss$thresholds[, c("Estimate", "Std. Error")], thr)
    expect_false("thresholds" %in% names(ss$coefficients))
    expect_output(print(ss), "Threshold coefficients:")
    expect_output(print(fit_ord), "Low\\|Medium = .*Medium\\|High = ")
    fit_pois <- glmmTMB(count ~ mined, family = poisson, data = Salamanders)
    expect_null(summary(fit_pois)$thresholds)
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
    ## silence warning about integer codes
    fit2 <- suppressWarnings(update(fit, data = dd2))
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
    expect_no_warning(glmmTMB(Sat ~ Infl, weights = Freq, data = housing,
                              family = ordinal()))
    ## a real REML fit still integrates out the fixed effects
    expect_equal(unique(names(fit$sdr$par.random)), "beta")
})

test_that("ordinal REML null model stays usable", {
    ## the ordinal intercept is always mapped, so y ~ 1 leaves no free
    ## fixed effect to integrate out and sdreport() returns no joint
    ## precision matrix; vcov() must fall back rather than error
    fit0 <- suppressWarnings(
        glmmTMB(Sat ~ 1, weights = Freq, data = housing,
                family = ordinal(), REML = TRUE))
    expect_length(fit0$sdr$par.random, 0)
    expect_null(fit0$sdr$jointPrecision)
    expect_s3_class(summary(fit0), "summary.glmmTMB")
    expect_true(is.matrix(vcov(fit0)$cond))
    expect_true(all(is.finite(confint(fit0)[, "Estimate"])))
    ## same for any family with no fixed effects at all (no 'cond' block to
    ## report there, but it must not error)
    dd <- data.frame(y = rnorm(50))
    expect_s3_class(vcov(glmmTMB(y ~ 0, data = dd, REML = TRUE)),
                    "vcov.glmmTMB")
})

## emmeans modes for ordinal fits (cumulative-link modes as for
## ordinal::clm / MASS::polr). Oracle records for the numeric assertions
## in the four blocks below (ID: type; asserting block/expectation; source):
##   ORD-EMM-1: live; "ordinal emmeans modes match ordinal::clm",
##     estimate + SE columns per (formula, mode); emmeans on ordinal::clm
##     (emmeans:::emm_basis.clm, ordinal package), logit and probit links
##   ORD-EMM-2: live; "ordinal emmeans on a mixed model", (a) prob column;
##     emmeans on ordinal::clmm (emmeans:::emm_basis.clmm)
##   ORD-EMM-3: invariant; same block (b) and "ordinal emmeans with a
##     mapped coefficient"; predict(., type = "probs") (an independent
##     route through the TMB template) must agree with the emmeans grid
##   ORD-EMM-4: closed-form; same block (c) P(Y <= j) = plogis(theta_j - eta)
##     and (d) E[class] = sum_j j * P(Y = j), written out in the test;
##     Agresti (2010) Analysis of Ordinal Categorical Data, ch. 3
test_that("ordinal emmeans modes match ordinal::clm", {
    skip_if_not_installed("emmeans")
    skip_if_not_installed("ordinal")
    data("wine", package = "ordinal")
    m_tmb <- glmmTMB(rating ~ temp + contact, family = ordinal, data = wine)
    m_clm <- ordinal::clm(rating ~ temp + contact, data = wine)
    m_tmb_p <- glmmTMB(rating ~ temp + contact,
                       family = ordinal(link = "probit"), data = wine)
    m_clm_p <- ordinal::clm(rating ~ temp + contact, data = wine,
                            link = "probit")
    cmp_clm <- function(fit_tmb, fit_clm, formula, mode) {
        s <- summary(emmeans::emmeans(fit_tmb, formula, mode = mode))
        s_clm <- summary(emmeans::emmeans(fit_clm, formula, mode = mode))
        en <- attr(s, "estName")
        expect_identical(en, attr(s_clm, "estName"), info = mode)
        expect_equal(s[[en]], s_clm[[en]], tolerance = 1e-4, info = mode)
        expect_equal(s[["SE"]], s_clm[["SE"]], tolerance = 1e-4, info = mode)
    }
    cases <- list(list(~ temp, "latent"),
                  list(~ cut | temp, "linear.predictor"),
                  list(~ cut | temp, "cum.prob"),
                  list(~ cut | temp, "exc.prob"),
                  list(~ rating | temp, "prob"),
                  list(~ temp, "mean.class"))
    for (cs in cases) cmp_clm(m_tmb, m_clm, cs[[1]], cs[[2]])
    cmp_clm(m_tmb_p, m_clm_p, ~ cut | temp, "cum.prob")
    ## default mode is latent, on which type = "response" is a no-op
    s_def <- summary(emmeans::emmeans(m_tmb, ~ temp))
    s_lat <- summary(emmeans::emmeans(m_tmb, ~ temp, mode = "latent"))
    expect_equal(s_def$emmean, s_lat$emmean)
    expect_equal(s_def$SE, s_lat$SE)
    s_resp <- summary(emmeans::emmeans(m_tmb, ~ temp, type = "response"))
    expect_equal(s_resp$emmean, s_def$emmean)
})

test_that("ordinal emmeans on a mixed model", {
    skip_if_not_installed("emmeans")
    skip_if_not_installed("ordinal")
    data("wine", package = "ordinal")
    m_mix <- glmmTMB(rating ~ temp + contact + (1 | judge),
                     family = ordinal, data = wine)
    m_clmm <- ordinal::clmm(rating ~ temp + contact + (1 | judge),
                            data = wine)
    ## (a) class probabilities agree with emmeans on ordinal::clmm
    s_prob <- summary(emmeans::emmeans(m_mix, ~ rating | temp + contact,
                                       mode = "prob"))
    s_clmm <- summary(emmeans::emmeans(m_clmm, ~ rating | temp + contact,
                                       mode = "prob"))
    expect_equal(s_prob$prob, s_clmm$prob, tolerance = 1e-3)
    s_cum <- summary(emmeans::emmeans(m_mix, ~ cut | temp + contact,
                                      mode = "cum.prob"))
    en_cum <- attr(s_cum, "estName")
    s_mc <- summary(emmeans::emmeans(m_mix, ~ temp + contact,
                                     mode = "mean.class"))
    theta <- family_params(m_mix)
    cells <- expand.grid(temp = levels(wine$temp),
                         contact = levels(wine$contact))
    for (i in seq_len(nrow(cells))) {
        cell <- cells[i, ]
        sel <- function(s) s$temp == cell$temp & s$contact == cell$contact
        ## (b) population-level predict() gives the same probabilities
        p <- drop(predict(m_mix, newdata = cell, type = "probs",
                          re.form = NA))
        expect_equal(s_prob$prob[sel(s_prob)], unname(p), tolerance = 1e-8)
        ## (c) cumulative probabilities: P(Y <= j) = plogis(theta_j - eta),
        ## eta = x'beta with the (always zero) intercept column removed
        Xc <- model.matrix(~ temp + contact, cell)
        Xc[, "(Intercept)"] <- 0
        eta <- drop(Xc %*% fixef(m_mix)$cond)
        expect_equal(s_cum[[en_cum]][sel(s_cum)],
                     unname(plogis(theta - eta)), tolerance = 1e-8)
        ## (d) mean class: E[Y] = sum_j j * P(Y = j)
        expect_equal(s_mc$mean.class[sel(s_mc)],
                     sum(seq_len(5) * p), tolerance = 1e-8)
    }
})

test_that("ordinal emmeans with a mapped coefficient", {
    skip_if_not_installed("emmeans")
    skip_if_not_installed("ordinal")
    data("wine", package = "ordinal")
    m_map <- glmmTMB(rating ~ temp + contact, family = ordinal, data = wine,
                     map = list(beta = factor(c(NA, 1, NA))),
                     start = list(beta = c(0, 0, 0)))
    expect_no_error(
        expect_no_warning(
            em <- emmeans::emmeans(m_map, ~ rating | temp + contact,
                                   mode = "prob")))
    s <- summary(em)
    cells <- expand.grid(temp = levels(wine$temp),
                         contact = levels(wine$contact))
    for (i in seq_len(nrow(cells))) {
        cell <- cells[i, ]
        p <- drop(predict(m_map, newdata = cell, type = "probs"))
        expect_equal(s$prob[s$temp == cell$temp & s$contact == cell$contact],
                     unname(p), tolerance = 1e-8)
    }
})

test_that("ordinal emmeans forces asymptotic ddf", {
    skip_if_not_installed("emmeans")
    skip_if_not_installed("ordinal")
    data("wine", package = "ordinal")
    m_mix <- glmmTMB(rating ~ temp + contact + (1 | judge),
                     family = ordinal, data = wine)
    msgs <- character(0)
    em <- withCallingHandlers(
        emmeans::emmeans(m_mix, ~ temp, ddf = "satterthwaite"),
        warning = function(w) {
            msgs <<- c(msgs, conditionMessage(w))
            invokeRestart("muffleWarning")
        })
    expect_length(msgs, 1L)
    expect_match(msgs, "using ddf 'asymptotic'")
    expect_true(all(summary(em)$df == Inf))
})
