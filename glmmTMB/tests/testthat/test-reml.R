stopifnot(require("testthat"),
          require("glmmTMB"),
          require("lme4"))


fm1.glmmTMB <- glmmTMB(Reaction ~ Days + (Days | Subject),
                       sleepstudy, REML=TRUE)

test_that("REML check against lmer", {
    ## Example 1: Compare results with lmer
    fm1.lmer <- lmer(Reaction ~ Days + (Days | Subject),
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
                 as.matrix(vcov(fm2.lmer)) , tolerance=1e-2)
})

test_that("REML with all parameters fixed", {
    john.alpha <- readRDS(system.file("test_data", "agridat_john.alpha.rds", package = "glmmTMB"))
    mod6_REML <- glmmTMB(yield ~ rep + (1 | gen),
                     start = list(theta =   log(sqrt(3)), betadisp =       log(5)),
                     map   = list(theta =  factor(c(NA)), betadisp = factor(c(NA))),
                     REML = TRUE,
                     data = john.alpha)
    mod6_ML <- update(mod6_REML, REML = FALSE)
    expect_equal(vcov(mod6_REML), vcov(mod6_ML))
})

test_that("REML with no free fixed effects stays usable", {
    ## with nothing random left at all -- no random effects, and no free
    ## fixed effects to integrate out -- sdreport() returns no
    ## jointPrecision, and vcov() used to fail with "attempt to set an
    ## attribute on NULL"
    set.seed(101)
    dd <- data.frame(y = rnorm(50), x = rnorm(50))
    m0 <- glmmTMB(y ~ 0, data = dd, REML = TRUE)
    expect_null(m0$sdr$jointPrecision)
    expect_s3_class(vcov(m0), "vcov.glmmTMB")
    expect_s3_class(summary(m0), "summary.glmmTMB")
    ## same when every fixed effect is map-fixed
    m1 <- glmmTMB(y ~ x, data = dd, REML = TRUE,
                  start = list(beta = c(0, 0)),
                  map = list(beta = factor(c(NA, NA))))
    expect_null(m1$sdr$jointPrecision)
    expect_s3_class(vcov(m1), "vcov.glmmTMB")
    ## a genuine REML fit is unaffected: beta really is integrated out
    expect_setequal(unique(names(fm1.glmmTMB$sdr$par.random)),
                    c("beta", "b"))
})

test_that("correct df.residual for REML=TRUE", {
    ## nobs = 180 - 6 (beta = 2 + betadisp = 1 + theta = 3)
    expect_equal(df.residual(fm1.glmmTMB), 174)
})
