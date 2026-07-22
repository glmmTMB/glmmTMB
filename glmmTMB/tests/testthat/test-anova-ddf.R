library(glmmTMB)

## compare anova.glmmTMB's joint Kenward-Roger/Satterthwaite F-tests
## (ddf != "asymptotic") against lme4 + pbkrtest, for (1) a multi-parameter
## (2 df) drop and (2) a chain of 3 nested models (2 sequential 1 df drops)
##
## NB: pbkrtest::SATmodcomp/KRmodcomp re-evaluate the lmer call internally
## (e.g. via get_devfun()/update()), so all such calls are made here at
## top level (where `dd` is visible for that re-evaluation) rather than
## inside test_that() blocks
if (requireNamespace("pbkrtest") && requireNamespace("lme4")) {

    set.seed(101)
    dd <- data.frame(x = rnorm(100), y = rnorm(100), f = factor(rep(1:10, each = 10)))
    dd$z <- simulate_new(~ x + y + (1|f), newdata = dd,
                          newparams = list(theta = 1, betadisp = -1, beta = c(0, 1, 1)),
                          family = gaussian)[[1]]
    ## pbkrtest::SATmodcomp/KRmodcomp re-fit internally via calls that get
    ## eval'd several frames up the call stack; putting `dd` in the global
    ## environment ensures it's found regardless of how this file is sourced
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

    rm(dd, envir = globalenv())
}
