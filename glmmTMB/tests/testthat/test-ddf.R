library(glmmTMB)
if (requireNamespace("pbkrtest") && requireNamespace("lme4")) {
    fm1_lmer <- lme4::lmer(formula(fm1), lme4::sleepstudy)
    fm2_lmer <- lme4::lmer(formula(fm2), lme4::sleepstudy)

    fm1 <- update(fm1, REML = TRUE)
    fm2 <- update(fm2, REML = TRUE)
    
    pbkrtest_dof <- function(m) {
        vva <- pbkrtest::vcovAdj(m)
        vv0 <- vcov(m)
        p <- length(fixef(m))
        Lmat <- diag(p)
        ddf <- apply(Lmat, 1, \(L) pbkrtest::ddf_Lb(vva, L, vv0))
        return(ddf)
    }

    test_that("approx KR df match with pbkrtest (fm1)", {
        df1 <- pbkrtest_dof(fm1_lmer)
        df2 <- dof_KR(fm1)
        expect_equal(df1, unname(c(df2)), tolerance = 1e-6)
    })

    test_that("approx KR df match with pbkrtest (fm2)", {
        df1 <- pbkrtest_dof(fm2_lmer)
        df2 <- dof_KR(fm2)
        expect_equal(df1, unname(c(df2)), tolerance = 1e-6)
    })

    test_that("KR in summary", {
        expect_identical(
            coef(summary(fm1, ddf = "kenward-roger"))$cond[,"ddf"],
            c(dof_KR(fm1))
        )
    })
}


test_that("Satt in summary", {
        expect_identical(
            unname(coef(summary(fm1, ddf = "satterthwaite"))$cond[,"ddf"]),
            c(dof_satt(fm1))
        )
    })               

test_that("Satterthwaite ddf match lmerTest (hard-coded reference values)", {
    ## reference values computed once via lmerTest 3.x (not at test time, to
    ## avoid adding a test dependency on lmerTest):
    ##   library(lmerTest)
    ##   fm1 <- lmerTest::lmer(Reaction ~ Days + (1 | Subject), sleepstudy, REML = TRUE)
    ##   fm2 <- lmerTest::lmer(Reaction ~ Days + (Days | Subject), sleepstudy, REML = TRUE)
    ##   coef(summary(fm1))[, "df"]; coef(summary(fm2))[, "df"]
    lmerTest_ddf_fm1 <- c(`(Intercept)` = 22.8102, Days = 161.0000)
    lmerTest_ddf_fm2 <- c(`(Intercept)` = 16.99973, Days = 16.99998)

    ## refit with REML = TRUE explicitly rather than relying on the
    ## conditional pbkrtest/lme4 block above (which only reassigns fm1/fm2 to
    ## REML fits if pbkrtest happens to be installed)
    fm1_reml <- update(fm1, REML = TRUE)
    fm2_reml <- update(fm2, REML = TRUE)

    ## tolerance of 1% comfortably covers the small, expected numerical
    ## difference between glmmTMB's finite-difference Satterthwaite
    ## implementation and lmerTest's (differences were ~0.03% in testing)
    expect_equal(unname(dof_satt(fm1_reml)), unname(lmerTest_ddf_fm1), tolerance = 1e-2)
    expect_equal(unname(dof_satt(fm2_reml)), unname(lmerTest_ddf_fm2), tolerance = 1e-2)
})

## emmeans

if (requireNamespace("emmeans")) {
    salamander1 <- up2date(readRDS(system.file("example_files", "salamander1.rds",
                                                package = "glmmTMB")))

    test_that("emmeans works with ddf='satterthwaite' (GH #1304)", {
        emm <- expect_no_error(
            suppressWarnings(emmeans::emmeans(salamander1, ~ mined, ddf = "satterthwaite"))
        )
        expect_true(all(is.finite(summary(emm)$df)))
    })

    test_that("satterthwaite dffun survives environment-stripping and per-contrast vector calls (GH #1304)", {
        rg <- suppressWarnings(emmeans::ref_grid(salamander1, ddf = "satterthwaite"))

        ## emmeans::ref_grid() replaces dffun's enclosing environment with
        ## baseenv() (glmmTMB#1304); rather than relying on that emmeans
        ## implementation detail, strip it here ourselves so the test keeps
        ## working even if emmeans stops doing so. dffun must not rely on
        ## any free variables (e.g. a captured copy of dof_satt) -- everything
        ## it needs must be reachable via the 'dfargs' argument instead of
        ## lexical scoping
        environment(rg@dffun) <- baseenv()

        ## dffun is called once per contrast with a bare vector (not a
        ## contrast matrix), so the wrapper must turn 'k' into a 1-row
        ## matrix before passing it on to dof_satt()
        k <- rg@linfct[1, ]
        expect_false(is.matrix(k))
        df_val <- suppressMessages(rg@dffun(k, rg@dfargs))
        expect_true(is.finite(df_val))
        expect_length(df_val, 1L)
    })

    ## constructed explicitly via update() (rather than relying on fm1's
    ## REML status, which is mutated earlier in this file if pbkrtest/lme4
    ## are available); update() reuses fm1's existing formula/data, so this
    ## doesn't need direct access to sleepstudy or the lme4 namespace
    fm1_ml_explicit <- update(fm1, REML = FALSE)
    fm1_reml_explicit <- update(fm1_ml_explicit, REML = TRUE)

    test_that("emmeans errors (rather than silently downgrading) for ddf='kenward-roger' on an ML fit", {
        expect_error(
            emmeans::emmeans(fm1_ml_explicit, "Days", ddf = "kenward-roger"),
            "requires a REML fit"
        )
    })

    test_that("emmeans allows ddf='kenward-roger' for a REML fit, with no warning", {
        emm <- expect_no_warning(emmeans::emmeans(fm1_reml_explicit, "Days", ddf = "kenward-roger"))
        expect_true(is.finite(summary(emm)$df[1]))
    })

    test_that("emmeans warns (does not error) for K-R/Satterthwaite on a non-Gaussian family", {
        expect_warning(
            emmeans::emmeans(salamander1, ~ mined, ddf = "satterthwaite"),
            "poorly understood"
        )
    })

    ## salamander1's family (poisson) has no estimated dispersion parameter
    ## (usesDispersion() == FALSE); check_ddf() used to hard-error for this
    ## unconditionally (contradicting the "warn, don't error" GLMM policy
    ## documented above and enforced by emmeans), so summary()/anova()/Anova()
    ## and emmeans() disagreed for exactly this case. satterthwaite's
    ## implementation doesn't need a dispersion parameter and should now warn
    ## (not error) consistently everywhere; kenward-roger's variance-component
    ## machinery genuinely doesn't support such families and should still
    ## error, but with a specific message rather than a crash
    test_that("summary() and emmeans() agree: ddf='satterthwaite' warns (does not error) for a family with no dispersion parameter", {
        expect_warning(summary(salamander1, ddf = "satterthwaite"), "poorly understood")
        expect_warning(emmeans::emmeans(salamander1, ~ mined, ddf = "satterthwaite"), "poorly understood")
    })

    salamander1_reml <- update(salamander1, REML = TRUE)

    test_that("summary() and emmeans() agree: ddf='kenward-roger' errors clearly (not an opaque eigen()/forceSymmetric crash) for a family with no dispersion parameter", {
        expect_error(summary(salamander1_reml, ddf = "kenward-roger"), "no estimated dispersion parameter")
        expect_error(emmeans::emmeans(salamander1_reml, ~ mined, ddf = "kenward-roger"), "no estimated dispersion parameter")
    })
}
