stopifnot(require("testthat"),
          require("glmmTMB"),
          require("lme4"))

data(sleepstudy, cbpp,
     package = "lme4")

context("variance structures")

## two equivalent diagonal constructions
fm1 <- glmmTMB(Reaction ~ Days + diag(Days| Subject), sleepstudy)
fm2 <- glmmTMB(Reaction ~ Days + ( 1  | Subject) + (0+Days | Subject),
               sleepstudy)
fm2L <- lmer(Reaction ~ Days + ( 1  | Subject) + (0+Days | Subject),
               sleepstudy, REML=FALSE)

fm3 <- glmmTMB(Reaction ~ Days + (Days| Subject), sleepstudy)
fm4 <- glmmTMB(Reaction ~ Days + cs(Days| Subject), sleepstudy)
fm3L <- lmer(Reaction ~ Days + ( Days  | Subject),
               sleepstudy, REML=FALSE)

test_that("diag", {
   ## two formulations of diag and lme4 all give same log-lik
   expect_equal(logLik(fm1),logLik(fm2L))
   expect_equal(logLik(fm1),logLik(fm2))
})

test_that("cs_us", {
## for a two-level factor, compound symmetry and unstructured
##  give same result
expect_equal(logLik(fm3),logLik(fm4))
expect_equal(logLik(fm3),logLik(fm3L))
})
