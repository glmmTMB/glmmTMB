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
##   ORD-EMM-4: closed-form; same block (c) P(Y <= j) = plogis(theta_j - eta),
##     written out in the test; Agresti (2010) Analysis of Ordinal
##     Categorical Data, ch. 3. (d) E[class] = sum_j j * P(Y = j) is the
##     definition of mean.class in vignette("models", package = "emmeans")
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

test_that("ordinal emmeans handles a rank-deficient fit", {
    skip_if_not_installed("emmeans")
    skip_if_not_installed("ordinal")
    data("wine", package = "ordinal")
    wine$dup <- wine$temp
    m_rd <- suppressWarnings(glmmTMB(rating ~ temp + contact + dup,
                                     family = ordinal, data = wine))
    expect_true(is.na(fixef(m_rd)$cond[["dupwarm"]]))
    em <- expect_no_error(
        emmeans::emmeans(m_rd, ~ rating | temp + contact, mode = "prob"))
    s <- summary(em)
    expect_false(anyNA(s$SE))
    cells <- expand.grid(temp = levels(wine$temp),
                         contact = levels(wine$contact))
    for (i in seq_len(nrow(cells))) {
        cell <- cells[i, ]
        cell$dup <- cell$temp
        p <- drop(predict(m_rd, newdata = cell, type = "probs"))
        expect_equal(s$prob[s$temp == cell$temp & s$contact == cell$contact],
                     unname(p), tolerance = 1e-8)
    }
})

test_that("ordinal emmeans carries a free intercept from a user map", {
    skip_if_not_installed("emmeans")
    skip_if_not_installed("ordinal")
    data("wine", package = "ordinal")
    ## a user map that leaves the intercept free (the default map fixes
    ## it to zero); the basis must include its estimate
    m_free <- glmmTMB(rating ~ temp + contact, family = ordinal, data = wine,
                      map = list(beta = factor(c(1, 2, 3))))
    expect_false(fixef(m_free)$cond[["(Intercept)"]] == 0)
    s <- summary(emmeans::emmeans(m_free, ~ rating | temp + contact,
                                  mode = "prob"))
    cells <- expand.grid(temp = levels(wine$temp),
                         contact = levels(wine$contact))
    for (i in seq_len(nrow(cells))) {
        cell <- cells[i, ]
        p <- drop(predict(m_free, newdata = cell, type = "probs"))
        expect_equal(s$prob[s$temp == cell$temp & s$contact == cell$contact],
                     unname(p), tolerance = 1e-8)
    }
})

test_that("ordinal emmeans latent mode honours rescale as MASS::polr", {
    skip_if_not_installed("emmeans")
    skip_if_not_installed("ordinal")
    skip_if_not_installed("MASS")
    data("wine", package = "ordinal")
    m_tmb <- glmmTMB(rating ~ temp + contact, family = ordinal, data = wine)
    m_polr <- MASS::polr(rating ~ temp + contact, data = wine, Hess = TRUE)
    s <- summary(emmeans::emmeans(m_tmb, ~ temp, mode = "latent",
                                  rescale = c(1, 10)))
    s_polr <- summary(emmeans::emmeans(m_polr, ~ temp, mode = "latent",
                                       rescale = c(1, 10)))
    expect_equal(s$emmean, s_polr$emmean, tolerance = 1e-4)
    expect_equal(s$SE, s_polr$SE, tolerance = 1e-4)
    s1 <- summary(emmeans::emmeans(m_tmb, ~ temp, mode = "latent"))
    expect_equal(s$emmean, 1 + 10 * s1$emmean, tolerance = 1e-8)
    expect_equal(s$SE, 10 * s1$SE, tolerance = 1e-8)
})

test_that("emmeans mode argument is rejected off the ordinal branch", {
    skip_if_not_installed("emmeans")
    skip_if_not_installed("ordinal")
    data("wine", package = "ordinal")
    m_num <- glmmTMB(as.numeric(rating) ~ temp + contact, data = wine)
    expect_error(emmeans::emmeans(m_num, ~ temp, mode = "prob"),
                 "only available for ordinal fits")
    m_ord <- glmmTMB(rating ~ temp + contact, family = ordinal, data = wine)
    expect_error(emmeans::emmeans(m_ord, ~ temp, mode = "nonsense"),
                 "'arg' should be one of")
})

## identifiability and reporting under a user map, an intercept-free
## formula, mapped psi and the profile-type interval methods. Oracle
## records for the numeric assertions (ID: type; asserting block; source):
##   ORD-ID-1: invariant; "ordinal user map on another component keeps the
##     intercept fixed" and "ordinal intercept-free formula is refitted
##     with an intercept"; the unmapped with-intercept fit fit_ord (GH #1348)
##   ORD-ID-2: live; "ordinal intercept-free formula is refitted with an
##     intercept", logLik and predict(type = "probs"); ordinal::clm on the
##     same intercept-free formula, which also assumes an intercept
##   ORD-ID-3: closed-form; "ordinal mapped psi: vcov, summary and Wald
##     confint", thresholds are a deterministic function of psi, so fixing
##     every psi element gives standard error 0 and an interval of width 0
##   ORD-ID-4: invariant; "ordinal profile-type intervals are labelled
##     psi", the uniroot Estimate column must equal the psi values stored
##     in the fitted object (two routes to the same number)
##   ORD-ID-5: invariant; "ordinal mapped psi: vcov, summary and Wald
##     confint", emmeans(mode = "prob") against predict(type = "probs")
##     on the same fit (two routes to the same probabilities)
##   ORD-ID-6: invariant; "ordinal profile-type intervals are labelled
##     psi", the profile interval against the uniroot interval at 1e-2
##     (two root finders on the same profile likelihood)
##   ORD-ID-7: invariant; "ordinal psi labels do not collide with a
##     fixed-effect column", the fixed-effect block of vcov(full = TRUE)
##     against vcov()$cond (two routes to the same matrix)
test_that("ordinal user map on another component keeps the intercept fixed", {
    ## a map on a different parameter vector must not disable the
    ## internal intercept map (partial matching of map$beta, GH #1348)
    fit_map <- glmmTMB(Sat ~ Infl + Type + Cont, weights = Freq,
                       data = housing, family = ordinal(),
                       map = list(betazi = factor()))
    expect_equal(fixef(fit_map)$cond[["(Intercept)"]], 0)
    bmap <- fit_map$obj$env$map[["beta"]]
    expect_true(is.na(bmap[[1]]))
    expect_false("(Intercept)" %in%
                 rownames(summary(fit_map)$coefficients$cond))
    expect_equal(fixef(fit_map)$cond, fixef(fit_ord)$cond, tolerance = 1e-6)
    expect_equal(family_params(fit_map), family_params(fit_ord),
                 tolerance = 1e-6)
    expect_equal(c(logLik(fit_map)), c(logLik(fit_ord)), tolerance = 1e-6)
})

test_that("ordinal intercept-free formula is refitted with an intercept", {
    ## as ordinal::clm and MASS::polr: warn, then assume the intercept
    for (ff in list(Sat ~ 0 + Infl + Type + Cont,
                    Sat ~ Infl + Type + Cont - 1)) {
        expect_warning(fit0 <- glmmTMB(ff, weights = Freq, data = housing,
                                       family = ordinal()),
                       "intercept is needed and assumed")
        expect_equal(fixef(fit0)$cond, fixef(fit_ord)$cond, tolerance = 1e-6)
        expect_equal(family_params(fit0), family_params(fit_ord),
                     tolerance = 1e-6)
        expect_equal(c(logLik(fit0)), c(logLik(fit_ord)), tolerance = 1e-6)
        ## the stored formula carries the intercept, so prediction on new
        ## data builds the same model matrix
        expect_equal(attr(terms(formula(fit0, component = "cond")),
                          "intercept"), 1L)
        expect_equal(predict(fit0, newdata = housing[1:6, ], type = "probs"),
                     predict(fit_ord, type = "probs")[1:6, ],
                     tolerance = 1e-6)
    }
    ## other families keep their intercept-free parameterization
    expect_no_warning(m <- glmmTMB(count ~ 0 + mined, data = Salamanders,
                                   family = poisson))
    expect_identical(names(fixef(m)$cond), c("minedyes", "minedno"))

    skip_if_not_installed("ordinal")
    ## live oracle: clm on the same intercept-free formula (it drops a
    ## different factor level, so compare the fit, not the coefficients)
    fit_clm <- suppressWarnings(
        ordinal::clm(Sat ~ 0 + Infl + Type + Cont, weights = Freq,
                     data = housing))
    fit0 <- suppressWarnings(glmmTMB(Sat ~ 0 + Infl + Type + Cont,
                                     weights = Freq, data = housing,
                                     family = ordinal()))
    expect_equal(c(logLik(fit0)), c(logLik(fit_clm)), tolerance = 1e-5)
    p_clm <- suppressWarnings(
        predict(fit_clm, newdata = housing[, c("Infl", "Type", "Cont")],
                type = "prob")$fit)
    expect_equal(unname(predict(fit0, type = "probs")), unname(p_clm),
                 tolerance = 1e-4)
    ## with a random effect
    data("wine", package = "ordinal")
    expect_warning(m0 <- glmmTMB(rating ~ 0 + temp + contact + (1 | judge),
                                 data = wine, family = ordinal()),
                   "intercept is needed")
    m1 <- glmmTMB(rating ~ temp + contact + (1 | judge), data = wine,
                  family = ordinal())
    expect_equal(c(logLik(m0)), c(logLik(m1)), tolerance = 1e-6)
    expect_equal(fixef(m0)$cond, fixef(m1)$cond, tolerance = 1e-5)
})

test_that("ordinal mapped psi: vcov, summary and Wald confint", {
    ## the thresholds are a joint function of all psi elements: fixing
    ## one of two leaves both thresholds free, fixing both makes them
    ## constants (here at psi = 0: equiprobable baseline categories)
    fit_p1 <- glmmTMB(Sat ~ Infl + Type + Cont, weights = Freq,
                      data = housing, family = ordinal(),
                      map = list(psi = factor(c(NA, 1))),
                      start = list(psi = c(0, 0)))
    fit_p2 <- glmmTMB(Sat ~ Infl + Type + Cont, weights = Freq,
                      data = housing, family = ordinal(),
                      map = list(psi = factor(c(NA, NA))),
                      start = list(psi = c(0, 0)))
    ## full vcov: NA rows and columns for the mapped psi entries only
    V1 <- vcov(fit_p1, full = TRUE)
    expect_true(all(c("psi1", "psi2") %in% rownames(V1)))
    expect_true(all(is.na(V1["psi1", ])) && all(is.na(V1[, "psi1"])))
    expect_true(is.finite(V1["psi2", "psi2"]))
    expect_true(is.finite(V1["InflHigh", "psi2"]))
    V2 <- vcov(fit_p2, full = TRUE)
    expect_true(all(is.na(V2[c("psi1", "psi2"), ])))
    expect_true(all(is.na(V2[, c("psi1", "psi2")])))
    expect_false(anyNA(V2[c("InflHigh", "ContHigh"), c("InflHigh", "ContHigh")]))
    for (fit in list(fit_p1, fit_p2)) {
        th <- summary(fit)$thresholds
        expect_true(all(is.finite(th[, c("Estimate", "Std. Error")])))
        ci <- confint(fit, component = "all")
        thr <- family_params(fit)
        expect_true(all(names(thr) %in% rownames(ci)))
        expect_true(all(is.finite(ci[names(thr), ])))
        ## the Wald half-width reproduces the summary standard error
        expect_equal(unname((ci[names(thr), 2] - ci[names(thr), 1]) /
                            (2 * qnorm(0.975))),
                     unname(th[, "Std. Error"]), tolerance = 1e-8)
    }
    expect_true(all(summary(fit_p1)$thresholds[, "Std. Error"] > 0))
    expect_equal(unname(summary(fit_p2)$thresholds[, "Std. Error"]), c(0, 0))
    ## psi = (0, 0) gives thresholds qlogis(1/3), qlogis(2/3)
    expect_equal(unname(family_params(fit_p2)), qlogis(c(1, 2) / 3))
    skip_if_not_installed("emmeans")
    s <- summary(emmeans::emmeans(fit_p2, ~ Sat | Infl + Type + Cont,
                                  mode = "prob"))
    expect_false(anyNA(s$SE))
    p <- drop(predict(fit_p2, newdata = data.frame(Infl = "Low", Type = "Tower",
                                                   Cont = "Low", Freq = 1),
                      type = "probs"))
    sel <- s$Infl == "Low" & s$Type == "Tower" & s$Cont == "Low"
    expect_equal(s$prob[sel], unname(p), tolerance = 1e-8)
})

test_that("ordinal profile-type intervals are labelled psi", {
    pars <- glmmTMB:::get_pars(fit_ord)
    psi_hat <- unname(pars[names(pars) == "psi"])
    ## the internal rows of the full vcov are psi1, psi2
    expect_identical(unname(tail(rownames(vcov(fit_ord, full = TRUE)), 2)),
                     c("psi1", "psi2"))
    ci_u <- confint(fit_ord, parm = "psi_", method = "uniroot")
    expect_identical(rownames(ci_u), c("psi1", "psi2"))
    ## the rows hold the psi values, not those of the parameters listed
    ## above them (the mapped intercept used to shift the labels)
    expect_equal(unname(ci_u[, "Estimate"]), psi_hat, tolerance = 1e-8)
    expect_true(all(ci_u[, 1] < psi_hat & psi_hat < ci_u[, 2]))
    pr <- profile(fit_ord, parm = "psi_", npts = 4)
    expect_identical(levels(pr$.par), c("psi1", "psi2"))
    ci_p <- confint(fit_ord, parm = "psi_", method = "profile", npts = 4)
    expect_identical(rownames(ci_p), c("psi1", "psi2"))
    expect_equal(unname(ci_p), unname(ci_u[, 1:2]), tolerance = 1e-2)
    ## Wald intervals stay on the threshold scale with threshold labels
    ci_w <- confint(fit_ord, component = "all")
    expect_true(all(c("Low|Medium", "Medium|High") %in% rownames(ci_w)))
    expect_false(any(grepl("^psi", rownames(ci_w))))
    expect_identical(rownames(summary(fit_ord)$thresholds),
                     c("Low|Medium", "Medium|High"))
})

test_that("ordinal intercept-free formula with a user beta map or start is an error", {
    ## the user's map or start vector is sized for the intercept-free
    ## model matrix; the added intercept would shift it
    expect_error(glmmTMB(Sat ~ 0 + Infl + Type + Cont, weights = Freq,
                         data = housing, family = ordinal(),
                         map = list(beta = factor(1:7))),
                 "needs an intercept in 'formula'")
    expect_error(glmmTMB(Sat ~ 0 + Infl + Type + Cont, weights = Freq,
                         data = housing, family = ordinal(),
                         start = list(beta = rep(0, 7))),
                 "needs an intercept in 'formula'")
    ## a one-sided formula (as predict() passes) keeps its terms when
    ## the intercept is added
    fit_old <- fit_ord
    fit_old$modelInfo$allForm$formula <- Sat ~ 0 + Infl + Type + Cont
    expect_warning(p_old <- predict(fit_old, newdata = housing[1:6, ],
                                    type = "probs"),
                   "intercept is needed")
    expect_equal(p_old, predict(fit_ord, type = "probs")[1:6, ],
                 tolerance = 1e-8)
})

test_that("ordinal threshold labels select the psi parameter in confint", {
    ci_w <- confint(fit_ord, component = "all")
    expect_equal(confint(fit_ord, parm = "Low|Medium"),
                 ci_w["Low|Medium", , drop = FALSE])
    expect_equal(confint(fit_ord, parm = c("ContHigh", "Medium|High")),
                 ci_w[c("cond.ContHigh", "Medium|High"), ],
                 ignore_attr = "dimnames")
    ci_u <- confint(fit_ord, parm = "psi_", method = "uniroot")
    expect_equal(confint(fit_ord, parm = "Medium|High", method = "uniroot"),
                 ci_u["psi2", , drop = FALSE])
})

test_that("ordinal psi labels do not collide with a fixed-effect column", {
    ## a factor named psi gives a fixed-effect column "psi2", the same
    ## label as the second internal threshold parameter; the full vcov
    ## must be filled by position, not by name
    skip_if_not_installed("ordinal")
    data("wine", package = "ordinal")
    set.seed(1)
    wine$psi <- factor(sample(1:2, nrow(wine), replace = TRUE))
    m <- glmmTMB(rating ~ temp + psi, data = wine, family = ordinal())
    V <- as.matrix(vcov(m, full = TRUE))
    expect_identical(unname(rownames(V)),
                     c("(Intercept)", "tempwarm", "psi2",
                       paste0("psi", 1:4)))
    ## every estimated row is finite; only the mapped intercept is NA
    expect_false(anyNA(V[-1, -1]))
    expect_true(all(is.na(V[1, ])))
    ## the fixed-effect block is the fixed-effect vcov
    expect_equal(unname(V[2:3, 2:3]),
                 unname(as.matrix(vcov(m, include_nonest = FALSE)$cond)))
    ## and the threshold standard errors come from the psi block
    th <- summary(m)$thresholds
    expect_true(all(is.finite(th[, "Std. Error"])))
    expect_true(all(th[, "Std. Error"] > 0))
})

test_that("ordinal working residuals are refused", {
    expect_error(residuals(fit_ord, type = "working"),
                 "working residuals are not defined for the ordinal family")
    ## the supported types still work
    expect_length(residuals(fit_ord, type = "response"), nrow(housing))
    expect_length(residuals(fit_ord, type = "dunn-smyth"), nrow(housing))
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
