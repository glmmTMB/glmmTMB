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
            suppressMessages(emmeans::emmeans(salamander1, ~ mined, ddf = "satterthwaite"))
        )
        expect_true(all(is.finite(summary(emm)$df)))
    })

    test_that("satterthwaite dffun survives environment-stripping and per-contrast vector calls (GH #1304)", {
        rg <- suppressMessages(emmeans::ref_grid(salamander1, ddf = "satterthwaite"))

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
}
