stopifnot(require("testthat"),
          require("glmmTMB"))

context("glmmTMBControl")

## Some selected L1-distances between two fits
distFits <- function(fit1, fit2) {
    s1 <- summary(fit1)
    s2 <- summary(fit2)
    glmmTMB:::namedList(
        max(abs((coef(s1)$cond - coef(s2)$cond)[,"Estimate"])),
        max(abs((coef(s1)$cond - coef(s2)$cond)[,"Std. Error"])),
        abs(logLik(fit1) - logLik(fit2))
    )
}

test_that("profile method", {

    myfit <- function(...) {
        glmmTMB(count ~ mined * spp + (1|site),
                family = poisson,
                data = Salamanders,
                control = glmmTMBControl(...))
    }

    m1 <- myfit( profile=FALSE )
    m2 <- myfit( profile=TRUE  )

    expect_true( all( distFits(m1, m2) < c(1e-4, 1e-2, 1e-4) ) )

    ## ###########################################################

    myfit <- function(...) {
        glmmTMB(count ~ mined * spp + (1|site),
                zi = ~ (1 | spp),
                family = poisson,
                data = Salamanders,
                control = glmmTMBControl(...))
    }

    m1 <- myfit( profile=FALSE )
    m2 <- myfit( profile=TRUE  )

    expect_true( all( distFits(m1, m2) < c(1e-4, 1e-2, 1e-4) ) )
})
