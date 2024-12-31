if (FALSE) {
    devtools::load_all("glmmTMB")
    devtools::document("glmmTMB")
    source("glmmTMB/tests/testthat/setup_makeex.R")
}
if (requireNamespace("pbkrtest") && requireNamespace("lme4")) {
    fm1_lmer <- lme4::lmer(formula(fm1), lme4::sleepstudy)
    fm2_lmer <- lme4::lmer(formula(fm2), lme4::sleepstudy)
    
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
        df2 <- dof_kenward(fm1)
        expect_equal(df1, unname(c(df2)), tolerance = 0.05)
    })

    test_that("approx KR df match with pbkrtest (fm2)", {
        df1 <- pbkrtest_dof(fm2_lmer)
        df2 <- dof_kenward(fm2)
        expect_equal(df1, unname(c(df2)), tolerance = 0.05)
    })

    dof_kenward(fm2NB)

}

