library(glmmTMB)

## compare anova.glmmTMB's joint Kenward-Roger/Satterthwaite F-tests
## (ddf != "asymptotic") against lme4 + pbkrtest, for (1) a multi-parameter
## (2 df) drop and (2) a chain of 3 nested models (2 sequential 1 df drops)
##
## NB: pbkrtest::SATmodcomp/KRmodcomp re-evaluate the lmer call internally
## (e.g. via get_devfun()/update()), so all such calls are made here at
## top level of `run_anova_ddf_tests()` (where `dd` is visible for that
## re-evaluation) rather than inside test_that() blocks. The whole thing is
## wrapped in a function (rather than left at file scope) so that on.exit()
## reliably restores the global environment even if an expectation below
## fails partway through -- on.exit() is a no-op at top level.
run_anova_ddf_tests <- function() {
    set.seed(101)
    dd <- data.frame(x = rnorm(100), y = rnorm(100), f = factor(rep(1:10, each = 10)))
    dd$z <- simulate_new(~ x + y + (1|f), newdata = dd,
                          newparams = list(theta = 1, betadisp = -1, beta = c(0, 1, 1)),
                          family = gaussian)[[1]]
    ## pbkrtest::SATmodcomp/KRmodcomp re-fit internally via calls that get
    ## eval'd several frames up the call stack; putting `dd` in the global
    ## environment ensures it's found regardless of how this file is sourced.
    ## Save/restore any pre-existing global `dd` so this is undone even if an
    ## expectation below fails partway through.
    had_dd <- exists("dd", envir = globalenv(), inherits = FALSE)
    if (had_dd) old_dd <- get("dd", envir = globalenv())
    on.exit({
        if (had_dd) {
            assign("dd", old_dd, envir = globalenv())
        } else if (exists("dd", envir = globalenv(), inherits = FALSE)) {
            rm("dd", envir = globalenv())
        }
    })
    assign("dd", dd, envir = globalenv())

    mnull <- glmmTMB(z ~ 1 + (1|f), data = dd, REML = TRUE)
    m0 <- glmmTMB(z ~ x + (1|f), data = dd, REML = TRUE)
    m1 <- glmmTMB(z ~ x + y + (1|f), data = dd, REML = TRUE)

    mnull_lmer <- lme4::lmer(z ~ 1 + (1|f), data = dd, REML = TRUE)
    m0_lmer <- lme4::lmer(z ~ x + (1|f), data = dd, REML = TRUE)
    m1_lmer <- lme4::lmer(z ~ x + y + (1|f), data = dd, REML = TRUE)

    ## multi-parameter (2 df) drop
    a_kr_drop <- anova(mnull, m1, ddf = "kenward-roger")
    kr_drop <- pbkrtest::KRmodcomp(m1_lmer, mnull_lmer)$test["Ftest", ]

    a_sat_drop <- anova(mnull, m1, ddf = "satterthwaite")
    sat_drop <- pbkrtest::SATmodcomp(m1_lmer, mnull_lmer)$test

    ## 3-model chain (2 sequential 1 df drops)
    a_kr_chain <- anova(mnull, m0, m1, ddf = "kenward-roger")
    kr_chain_1 <- pbkrtest::KRmodcomp(m0_lmer, mnull_lmer)$test["Ftest", ]
    kr_chain_2 <- pbkrtest::KRmodcomp(m1_lmer, m0_lmer)$test["Ftest", ]

    a_sat_chain <- anova(mnull, m0, m1, ddf = "satterthwaite")
    sat_chain_1 <- pbkrtest::SATmodcomp(m0_lmer, mnull_lmer)$test
    sat_chain_2 <- pbkrtest::SATmodcomp(m1_lmer, m0_lmer)$test

    ## small numerical differences are expected: glmmTMB and lme4 use
    ## different optimizers/parameterizations for the same REML fit
    ftol <- 1e-3

    test_that("anova KR F-test matches pbkrtest::KRmodcomp (multi-parameter drop)", {
        expect_equal(a_kr_drop$F[2], kr_drop[["stat"]], tolerance = ftol)
        expect_equal(a_kr_drop$`Num Df`[2], kr_drop[["ndf"]], tolerance = ftol)
        expect_equal(a_kr_drop$`Den Df`[2], kr_drop[["ddf"]], tolerance = ftol)
        expect_equal(a_kr_drop$`Pr(>F)`[2], kr_drop[["p.value"]], tolerance = ftol)
    })

    test_that("anova Satterthwaite F-test matches pbkrtest::SATmodcomp (multi-parameter drop)", {
        expect_equal(a_sat_drop$F[2], sat_drop$statistic, tolerance = ftol)
        expect_equal(a_sat_drop$`Num Df`[2], sat_drop$ndf, tolerance = ftol)
        expect_equal(a_sat_drop$`Den Df`[2], sat_drop$ddf, tolerance = ftol)
        expect_equal(a_sat_drop$`Pr(>F)`[2], sat_drop$p.value, tolerance = ftol)
    })

    test_that("anova KR F-test matches pbkrtest::KRmodcomp (3-model chain)", {
        expect_equal(a_kr_chain$F[2], kr_chain_1[["stat"]], tolerance = ftol)
        expect_equal(a_kr_chain$`Den Df`[2], kr_chain_1[["ddf"]], tolerance = ftol)
        expect_equal(a_kr_chain$`Pr(>F)`[2], kr_chain_1[["p.value"]], tolerance = ftol)
        expect_equal(a_kr_chain$F[3], kr_chain_2[["stat"]], tolerance = ftol)
        expect_equal(a_kr_chain$`Den Df`[3], kr_chain_2[["ddf"]], tolerance = ftol)
        expect_equal(a_kr_chain$`Pr(>F)`[3], kr_chain_2[["p.value"]], tolerance = ftol)
    })

    test_that("anova Satterthwaite F-test matches pbkrtest::SATmodcomp (3-model chain)", {
        expect_equal(a_sat_chain$F[2], sat_chain_1$statistic, tolerance = ftol)
        expect_equal(a_sat_chain$`Den Df`[2], sat_chain_1$ddf, tolerance = ftol)
        expect_equal(a_sat_chain$`Pr(>F)`[2], sat_chain_1$p.value, tolerance = ftol)
        expect_equal(a_sat_chain$F[3], sat_chain_2$statistic, tolerance = ftol)
        expect_equal(a_sat_chain$`Den Df`[3], sat_chain_2$ddf, tolerance = ftol)
        expect_equal(a_sat_chain$`Pr(>F)`[3], sat_chain_2$p.value, tolerance = ftol)
    })

    ## Anova() (car Type II/III) reuses the same joint K-R/Satterthwaite
    ## machinery via a hypothesis matrix rather than a pair of models; for an
    ## additive-only (no-interaction) model, Type II's per-term F-test for
    ## "x" is exactly the same hypothesis as dropping "x" while keeping "y",
    ## so it should agree with dropping x from m1 above via pbkrtest directly
    if (requireNamespace("car", quietly = TRUE)) {
        m_no_x <- glmmTMB(z ~ y + (1|f), data = dd, REML = TRUE)
        m_no_x_lmer <- lme4::lmer(z ~ y + (1|f), data = dd, REML = TRUE)

        kr_x <- pbkrtest::KRmodcomp(m1_lmer, m_no_x_lmer)$test["Ftest", ]
        sat_x <- pbkrtest::SATmodcomp(m1_lmer, m_no_x_lmer)$test

        Anova_kr <- car::Anova(m1, ddf = "kenward-roger")
        Anova_sat <- car::Anova(m1, ddf = "satterthwaite")
        Anova_sat_III <- car::Anova(m1, type = "III", ddf = "satterthwaite")

        test_that(".satt_adjust_joint() reuses the .satt_precompute() cache without changing results", {
            ## a fresh copy of m1 so this doesn't depend on cache state left
            ## behind by other tests -- the cache lives on model$obj$env, a
            ## genuine (shared-by-reference) environment, not a per-copy list
            m1_fresh <- glmmTMB(z ~ x + y + (1|f), data = dd, REML = TRUE)
            expect_null(m1_fresh$obj$env$.satt_cache)
            res_cold <- car::Anova(m1_fresh, ddf = "satterthwaite")
            expect_false(is.null(m1_fresh$obj$env$.satt_cache))
            res_warm <- car::Anova(m1_fresh, ddf = "satterthwaite")
            expect_identical(res_cold[["F"]], res_warm[["F"]])
            expect_identical(res_cold[["Den Df"]], res_warm[["Den Df"]])
        })

        test_that("Anova() Type II KR F-test matches pbkrtest::KRmodcomp (single-term drop)", {
            expect_equal(Anova_kr["x", "F"], kr_x[["stat"]], tolerance = ftol)
            expect_equal(Anova_kr["x", "Num Df"], kr_x[["ndf"]], tolerance = ftol)
            expect_equal(Anova_kr["x", "Den Df"], kr_x[["ddf"]], tolerance = ftol)
            expect_equal(Anova_kr["x", "Pr(>F)"], kr_x[["p.value"]], tolerance = ftol)
        })

        test_that("Anova() Type II Satterthwaite F-test matches pbkrtest::SATmodcomp (single-term drop)", {
            expect_equal(Anova_sat["x", "F"], sat_x$statistic, tolerance = ftol)
            expect_equal(Anova_sat["x", "Den Df"], sat_x$ddf, tolerance = ftol)
            expect_equal(Anova_sat["x", "Pr(>F)"], sat_x$p.value, tolerance = ftol)
        })

        test_that("Anova() Type II and Type III agree for an additive (no-interaction) model", {
            ## with no interactions and no other terms sharing marginality
            ## with "x", Type II and Type III should give the same test for "x"
            expect_equal(Anova_sat["x", "F"], Anova_sat_III["x", "F"], tolerance = 1e-8)
            expect_equal(Anova_sat["x", "Den Df"], Anova_sat_III["x", "Den Df"], tolerance = 1e-8)
        })
    }
}

if (requireNamespace("pbkrtest") && requireNamespace("lme4")) {
    run_anova_ddf_tests()
}

## behavioral checks that don't need pbkrtest/lme4 cross-validation: the
## "same checks" (REML requirement, GLMM warning, no-random-effects
## fallback, unsupported combinations) shared with summary()/anova()/emmeans()
if (requireNamespace("car", quietly = TRUE)) {
    set.seed(303)
    n_g <- 12
    g_sizes <- sample(5:25, n_g, replace = TRUE)
    dd_anova <- do.call(rbind, lapply(seq_len(n_g), function(gi) {
        n <- g_sizes[gi]
        data.frame(g = factor(gi), f = factor(sample(LETTERS[1:4], n, replace = TRUE)))
    }))
    dd_anova$y_gauss <- rnorm(nrow(dd_anova))
    dd_anova$y_nb <- rnbinom(nrow(dd_anova), mu = 5, size = 2)
    dd_anova$y_pois <- rpois(nrow(dd_anova), lambda = 5)

    m_ml <- glmmTMB(y_gauss ~ f + (1|g), data = dd_anova, family = gaussian, REML = FALSE)
    m_reml <- update(m_ml, REML = TRUE)
    m_nb <- glmmTMB(y_nb ~ f + (1|g), data = dd_anova, family = nbinom2, REML = TRUE)
    m_norandom <- glmmTMB(y_gauss ~ f, data = dd_anova, family = gaussian)
    ## poisson has no estimated dispersion parameter (usesDispersion() ==
    ## FALSE), unlike m_nb (nbinom2) above
    m_pois <- glmmTMB(y_pois ~ f + (1|g), data = dd_anova, family = poisson, REML = TRUE)

    test_that("Anova() ddf='kenward-roger' hard-errors (does not silently downgrade) on an ML fit", {
        expect_error(car::Anova(m_ml, ddf = "kenward-roger"), "requires a REML fit")
    })

    test_that("Anova() ddf='satterthwaite' works on an ML fit (no REML requirement)", {
        expect_no_error(expect_no_warning(car::Anova(m_ml, ddf = "satterthwaite")))
    })

    test_that("Anova() warns (does not error) for ddf on a non-Gaussian family", {
        expect_warning(car::Anova(m_nb, ddf = "satterthwaite"), "poorly understood")
        expect_warning(car::Anova(m_nb, ddf = "kenward-roger"), "poorly understood")
    })

    test_that("Anova() falls back to residual df (no error/warning about REML) when there are no random effects", {
        res <- expect_no_warning(suppressMessages(car::Anova(m_norandom, ddf = "kenward-roger")))
        expect_equal(unname(res["f", "Den Df"]), unname(df.residual(m_norandom)))
    })

    test_that("Anova() rejects ddf combined with a user-supplied vcov.", {
        expect_error(
            car::Anova(m_reml, ddf = "kenward-roger", vcov. = vcov(m_reml)$cond),
            "user-supplied 'vcov.'"
        )
    })

    test_that("Anova() rejects ddf for component != 'cond'", {
        ## convergence quality is irrelevant here -- the ddf/component guard
        ## fires before any numerical work that depends on it
        m_zi <- suppressWarnings(
            glmmTMB(y_nb ~ f + (1|g), ziformula = ~ f, data = dd_anova, family = nbinom2)
        )
        expect_error(car::Anova(m_zi, component = "zi", ddf = "satterthwaite"),
                    "only supported for component")
    })

    test_that("Anova() rejects ddf for models with map-fixed conditional coefficients", {
        m_map <- glmmTMB(y_gauss ~ f + (1|g), data = dd_anova, REML = TRUE,
                         map = list(beta = factor(c(1, NA, 2, 3))))
        expect_error(car::Anova(m_map, ddf = "kenward-roger"), "map-fixed")
    })

    test_that("Anova() test.statistic='F' without ddf is a clear error", {
        expect_error(car::Anova(m_reml, test.statistic = "F"), "ddf=")
    })

    test_that("Anova() rejects an explicit test.statistic='Chisq' combined with ddf != 'asymptotic'", {
        expect_error(
            car::Anova(m_reml, ddf = "kenward-roger", test.statistic = "Chisq"),
            "test.statistic='Chisq' cannot be combined with ddf"
        )
    })

    test_that("Anova() silently gives an F table when ddf is set without specifying test.statistic", {
        res <- expect_no_warning(expect_no_error(car::Anova(m_reml, ddf = "kenward-roger")))
        expect_true("F" %in% names(res))
        expect_match(attr(res, "heading")[1], "F tests")
    })

    test_that("Anova() ddf='kenward-roger' errors clearly (not an opaque eigen()/forceSymmetric crash) for a family with no dispersion parameter", {
        expect_error(
            car::Anova(m_pois, ddf = "kenward-roger"),
            "no estimated dispersion parameter"
        )
    })

    test_that("Anova() ddf='satterthwaite' works (warning, not error) for a family with no dispersion parameter", {
        res <- expect_warning(car::Anova(m_pois, ddf = "satterthwaite"), "poorly understood")
        expect_true(all(is.finite(res[["Den Df"]])))
    })
}
