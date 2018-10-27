stopifnot(require("testthat"),
          require("glmmTMB"),
          require("lme4"))

context("REML")

test_that("REML check against lmer", {
    ## Example 1: Compare results with lmer
    fm1.lmer <- lmer(Reaction ~ Days + (Days | Subject),
                     sleepstudy, REML=TRUE)
    fm1.glmmTMB <- glmmTMB(Reaction ~ Days + (Days | Subject),
                           sleepstudy, REML=TRUE)
    expect_equal( logLik(fm1.lmer) , logLik(fm1.glmmTMB) )
    expect_equal(as.vector(predict(fm1.lmer)) ,
                 predict(fm1.glmmTMB), tolerance=2e-3)
    expect_equal(vcov(fm1.glmmTMB)$cond,
                 as.matrix(vcov(fm1.lmer)) , tolerance=1e-3)
    ## Example 2: Compare results with lmer
    data(Orthodont,package="nlme")
    Orthodont$nsex <- as.numeric(Orthodont$Sex=="Male")
    Orthodont$nsexage <- with(Orthodont, nsex*age)
    fm2.lmer <- lmer(distance ~ age + (age|Subject) + (0+nsex|Subject) +
                         (0 + nsexage|Subject), data=Orthodont, REML=TRUE,
          control=lmerControl(check.conv.grad = .makeCC("warning", tol = 5e-3)))
    fm2.glmmTMB <- glmmTMB(distance ~ age + (age|Subject) + (0+nsex|Subject) +
                               (0 + nsexage|Subject), data=Orthodont, REML=TRUE)
    expect_equal( logLik(fm2.lmer) , logLik(fm2.glmmTMB), tolerance=1e-5 )
    expect_equal(as.vector(predict(fm2.lmer)) ,
                 predict(fm2.glmmTMB), tolerance=1e-4)
    expect_equal(vcov(fm2.glmmTMB)$cond,
                 as.matrix(vcov(fm2.lmer)) , tolerance=1e-3)
})
