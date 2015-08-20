stopifnot(require("testthat"),
          require("glmmTMB"))

data(sleepstudy, cbpp,
     package = "lme4")
cbpp <- transform(cbpp,prop=incidence/size)

context("Very basic glmmTMB fitting")


test_that("Basic Gaussian Sleepdata examples", {
    expect_is(fm0 <- glmmTMB(Reaction ~ 1    + ( 1  | Subject), sleepstudy), "glmmTMB")
    expect_is(fm1 <- glmmTMB(Reaction ~ Days + ( 1  | Subject), sleepstudy), "glmmTMB")
    expect_is(fm2 <- glmmTMB(Reaction ~ Days + (Days| Subject), sleepstudy), "glmmTMB")
    expect_is(fm3 <- glmmTMB(Reaction ~ Days + ( 1  | Subject) + (0+Days | Subject),
                             sleepstudy), "glmmTMB")

    ## *_equivalent(.), a misnomer, means  all.equal(*, check.attributes = FALSE):
    ## but can't use tolerance with *_equivalent (sigh) ...
    expect_equal(unname(fixef(fm0)), 298.508, tolerance = .0001)
    expect_equal(fixef(fm0), c("(Intercept)" = 298.508), tolerance = .0001)

    expect_equal(unname(fixef(fm1)), c(251.405, 10.4673),   tolerance = .0001)
    expect_equal(fixef(fm1), c("(Intercept)" = 251.405, Days = 10.4673),
                 tolerance = .0001)
})

test_that("Update Gaussian", {
  ## call doesn't match (formula gets mangled?)
  ## timing different
  ## FIXME: more redundancy
  fm0 <- glmmTMB(Reaction ~ 1    + ( 1  | Subject), sleepstudy)
  fm1 <- glmmTMB(Reaction ~ Days + ( 1  | Subject), sleepstudy)
  fm1u <- update(fm0, . ~ . + Days)
  fm1u$call <- fm1$call
  fm1u$optTime <- fm1$optTime
    expect_equal(fm1,fm1u )
})


test_that("Variance structures", {
  ## FIXME, redundant
  expect_is(fm2 <- glmmTMB(Reaction ~ Days + (Days| Subject), sleepstudy), 
            "glmmTMB")
  expect_is(fm2us <- glmmTMB(Reaction ~ Days + us(Days| Subject), sleepstudy), 
            "glmmTMB")
  expect_is(fm2cs <- glmmTMB(Reaction ~ Days + cs(Days| Subject), sleepstudy), 
            "glmmTMB")
  expect_is(fm2diag <- glmmTMB(Reaction ~ Days + diag(Days| Subject), sleepstudy), 
            "glmmTMB")
  gm <- glmmTMB:::getME
  expect_equal(gm(fm2,"theta"),gm(fm2us,"theta"))
  ## FIXME: more here, compare results against lme4 ...
})

test_that("Sleepdata Variance components", {
    ## TODO: Variance Components ("theta"s)
})


test_that("Basic Binomial CBPP examples", {
    ## intercept-only fixed effect
    expect_is(gm0 <- glmmTMB(prop ~ 1 + (1|herd),
                             weights=size,
                             data = cbpp, family=binomial()), "glmmTMB")
    expect_is(gm1 <- glmmTMB(prop ~ 1 + period + (1|herd),
                             weights=size,
                             data = cbpp, family=binomial()), "glmmTMB")

    expect_equal(fixef(gm0), c("(Intercept)" = -2.04567), tolerance = .0005)
    expect_equal(fixef(gm1), c("(Intercept)" = -1.39834,
                               period2 = -0.991925, period3 = -1.12822, 
                               period4 = -1.57975),
                  tolerance = .001) # <- TODO: lower eventually
})


test_that("Update Binomial", {
  ## call doesn't match (formula gets mangled?)
  ## timing different
  gm0 <- glmmTMB(prop ~ 1 + (1|herd),
                             weights=size,
                             data = cbpp, family=binomial())
  gm1 <- glmmTMB(prop ~ period + (1|herd),
                 weights=size,
                 data = cbpp, family=binomial())
  gm1u <- update(gm0, . ~ . + period)
  gm1u$call <- gm1$call
  gm1u$optTime <- gm1$optTime
    expect_equal(gm1, gm1u)
})




if(require("lme4")) {

    L <- load(system.file("testdata", "lme-tst-fits.rda",
                          package="lme4", mustWork=TRUE))
    message("loaded testdata from lme4:\n ",
            strwrap(paste(L, collapse = ", ")))

    if(FALSE) { ## part of the above [not recreated here for speed mostly:]
        ## intercept only in both fixed and random effects
        fit_sleepstudy_0 <- lmer(Reaction ~ 1 + (1|Subject), sleepstudy)
        ## fixed slope, intercept-only RE
        fit_sleepstudy_1 <- lmer(Reaction ~ Days + (1|Subject), sleepstudy)
        ## fixed slope, intercept & slope RE
        fit_sleepstudy_2 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
        ## fixed slope, independent intercept & slope RE
        fit_sleepstudy_3 <- lmer(Reaction ~ Days + (1|Subject)+ (0+Days|Subject), sleepstudy)

        cbpp$obs <- factor(seq(nrow(cbpp)))
        ## intercept-only fixed effect
        fit_cbpp_0 <- glmer(cbind(incidence, size-incidence) ~ 1 + (1|herd),
                            cbpp, family=binomial)
        ## include fixed effect of period
        fit_cbpp_1 <- update(fit_cbpp_0, . ~ . + period)
        ## include observation-level RE
        fit_cbpp_2 <- update(fit_cbpp_1, . ~ . + (1|obs))
        ## specify formula by proportion/weights instead
        fit_cbpp_3 <- update(fit_cbpp_1, incidence/size ~ period + (1 | herd), weights = size)
    }

    ## Now check closeness to lme4 results


}# comparing with lme4 if there..
