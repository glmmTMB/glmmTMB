stopifnot(require("testthat"),
          require("glmmTMB"))

data(sleepstudy, cbpp,
     package = "lme4")

data(quine, package="MASS")

## n.b. for test_that, this must be assigned within the global
## environment ...

cbpp <<- transform(cbpp, prop = incidence/size, obs=factor(seq(nrow(cbpp))))

## utility: hack/replace parts of the updated result that will
##  be cosmetically different
matchForm <- function(obj, objU, family=FALSE) {
  for(cmp in c("call","frame")) # <- more?
     objU[[cmp]] <- obj[[cmp]]
     ## Q: why are formulas equivalent but not identical?  A: their environments may differ
  objU$modelInfo$allForm <- obj$modelInfo$allForm
  if (family)  objU$modelInfo$family <- obj$modelInfo$family
  return(objU)
}

context("Very basic glmmTMB fitting")

lm0 <- lm(Reaction~Days,sleepstudy)
fm00 <- glmmTMB(Reaction ~ Days, sleepstudy)
fm0 <- glmmTMB(Reaction ~ 1    + ( 1  | Subject), sleepstudy)
fm1 <- glmmTMB(Reaction ~ Days + ( 1  | Subject), sleepstudy)
fm2 <- glmmTMB(Reaction ~ Days + (Days| Subject), sleepstudy)
fm3 <- glmmTMB(Reaction ~ Days + ( 1  | Subject) + (0+Days | Subject),
               sleepstudy)

test_that("Basic Gaussian Sleepdata examples", {
    expect_is(fm00, "glmmTMB")
    expect_is(fm0, "glmmTMB")
    expect_is(fm1, "glmmTMB")
    expect_is(fm2, "glmmTMB")
    expect_is(fm3, "glmmTMB")

    expect_equal(fixef(fm00)[[1]],coef(lm0),tol=1e-5)
    expect_equal(sigma(fm00)*sqrt(nobs(fm00)/(df.residual(fm00)+1)),
                 summary(lm0)$sigma,tol=1e-5)
    expect_equal(fixef(fm0)[[1]], c("(Intercept)" = 298.508), tolerance = .0001)
    expect_equal(fixef(fm1)[[1]], c("(Intercept)" = 251.405, Days = 10.4673),
                 tolerance = .0001)
    expect_equal(fixef(fm2)$cond, fixef(fm1)$cond, tolerance = 1e-5)# seen 1.042 e-6
    expect_equal(fixef(fm3)$cond, fixef(fm1)$cond, tolerance = 5e-6)# seen 2.250 e-7

    expect_equal(head(ranef(fm0)$cond$Subject[,1],3),
                 c(37.4881849228705, -71.5589277273216, -58.009085500647),
                 tolerance=1e-5)
    ## test *existence* of summary method -- nothing else for now
    expect_is(suppressWarnings(summary(fm3)),"summary.glmmTMB")
})

test_that("Update Gaussian", {
  ## call doesn't match (formula gets mangled?)
  ## timing different
  fm1u <- update(fm0, . ~ . + Days)
  expect_equal(fm1, matchForm(fm1, fm1u))
})


test_that("Variance structures", {
  ## above: fm2     <- glmmTMB(Reaction ~ Days +     (Days| Subject), sleepstudy)
  expect_is(fm2us   <- glmmTMB(Reaction ~ Days +   us(Days| Subject), sleepstudy), "glmmTMB")
  expect_is(fm2cs   <- glmmTMB(Reaction ~ Days +   cs(Days| Subject), sleepstudy), "glmmTMB")
  expect_is(fm2diag <- glmmTMB(Reaction ~ Days + diag(Days| Subject), sleepstudy), "glmmTMB")
  expect_equal(getME(fm2,  "theta"),
               getME(fm2us,"theta"))
  ## FIXME: more here, compare results against lme4 ...
})

test_that("Sleepdata Variance components", {
    expect_equal(c(unlist(VarCorr(fm3))),
                 c(cond.Subject = 584.247907378213, cond.Subject.1 = 33.6332741779585),
                 tolerance=1e-5)
})

gm0 <<- glmmTMB(cbind(incidence, size-incidence) ~ 1 +      (1|herd),
               data = cbpp, family=binomial())
gm1 <<- glmmTMB(cbind(incidence, size-incidence) ~ period + (1|herd),
               data = cbpp, family=binomial())

test_that("Basic Binomial CBPP examples", {

    ## Basic Binomial CBPP examples ---- intercept-only fixed effect
    expect_is(gm0, "glmmTMB")
    expect_is(gm1, "glmmTMB")
    expect_equal(fixef(gm0)[[1]], c("(Intercept)" = -2.045671), tolerance = 1e-3)#lme4 results
    expect_equal(fixef(gm1)[[1]], c("(Intercept)" = -1.398343,#lme4 results
                               period2 = -0.991925, period3 = -1.128216,
                               period4 = -1.579745),
                  tolerance = 1e-3) # <- TODO: lower eventually

})

test_that("Multiple RE, reordering", {
    ### Multiple RE,  reordering

    tmb1 <- glmmTMB(cbind(incidence, size-incidence) ~ period + (1|herd) + (1|obs),
                    data = cbpp, family=binomial())
    tmb2 <- glmmTMB(cbind(incidence, size-incidence) ~ period + (1|obs) + (1|herd),
                    data = cbpp, family=binomial())
    expect_equal(fixef(tmb1), fixef(tmb2),                   tolerance = 1e-8)
    expect_equal(getME(tmb1, "theta"), getME(tmb2, "theta")[c(2,1)], tolerance = 5e-7)
})

test_that("Alternative family specifications [via update(.)]", {
    ## intercept-only fixed effect

    res_chr <- matchForm(gm0, update(gm0, family= "binomial"))
    expect_equal(gm0, res_chr)
    expect_equal(gm0, matchForm(gm0, update(gm0, family= binomial())))
    expect_warning(res_list <- matchForm(gm0, update(gm0, family= list(family = "binomial",
                                                       link = "logit")),
                                         family=TRUE))
    expect_equal(gm0, res_list)
})

test_that("Update Binomial", {
  ## matchForm(): call doesn't match (formula gets mangled?)
  ## timing different
  gm1u <- update(gm0, . ~ . + period)
  expect_equal(gm1, matchForm(gm1, gm1u))
})

test_that("internal structures", {
  ## RE terms only in cond and zi model, not disp: GH #79
  expect_equal(names(fm0$modelInfo$reTrms),
               c("cond","zi"))
})

test_that("close to lme4 results", {
    expect_true(require("lme4"))
    L <- load(system.file("testdata", "lme-tst-fits.rda",
                          package="lme4", mustWork=TRUE))
    expect_is(L, "character")
    message("Loaded testdata from lme4:\n ",
            paste(strwrap(paste(L, collapse = ", ")),
                  collapse = "\n "))

    if(FALSE) { ## part of the above [not recreated here for speed mostly:]
        ## intercept only in both fixed and random effects
        fit_sleepstudy_0 <- lmer(Reaction ~   1  + ( 1 | Subject), sleepstudy)
        ## fixed slope, intercept-only RE
        fit_sleepstudy_1 <- lmer(Reaction ~ Days + ( 1 | Subject), sleepstudy)
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

    ## What we really want to compare against - Maximum Likelihood (package 'DESCRIPTION' !)
    fi_0 <- lmer(Reaction ~   1  + ( 1  | Subject), sleepstudy, REML=FALSE)
    fi_1 <- lmer(Reaction ~ Days + ( 1  | Subject), sleepstudy, REML=FALSE)
    fi_2 <- lmer(Reaction ~ Days + (Days| Subject), sleepstudy, REML=FALSE)
    fi_3 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject),
                 sleepstudy, REML=FALSE)

    ## Now check closeness to lme4 results

    ## ......................................
})

context("trickier examples")

data(Owls)
## is <<- necessary ... ?
Owls <- transform(Owls,
                   ArrivalTime=scale(ArrivalTime,center=TRUE,scale=FALSE),
                   NCalls= SiblingNegotiation) 

test_that("basic zero inflation", {
	expect_true(require("pscl"))
	o0.tmb <- glmmTMB(NCalls~(FoodTreatment + ArrivalTime) * SexParent + 
                              offset(logBroodSize), 
                          ziformula=~1, data = Owls,
                          family=poisson(link = "log"))
	o0.pscl <-zeroinfl(NCalls~(FoodTreatment + ArrivalTime) * SexParent + 
        offset(logBroodSize)|1, data = Owls)
    expect_equal(summary(o0.pscl)$coefficients$count, summary(o0.tmb)$coefficients$cond, tolerance=1e-5)
    expect_equal(summary(o0.pscl)$coefficients$zero, summary(o0.tmb)$coefficients$zi, tolerance=1e-5)
 
    o1.tmb <- glmmTMB(NCalls~(FoodTreatment + ArrivalTime) * SexParent + 
        offset(logBroodSize) + diag(1 | Nest), 
        ziformula=~1, data = Owls, family=poisson(link = "log"))
	expect_equal(ranef(o1.tmb)$cond$Nest[1,1], -0.484, tolerance=1e-2) #glmmADMB gave -0.4842771
})

test_that("alternative binomial model specifications", {
    d <<- data.frame(y=1:10,N=20,x=1) ## n.b. global assignment for testthat
    m0 <- suppressWarnings(glmmTMB(cbind(y,N-y) ~ 1, data=d, family=binomial()))
    m3 <- glmmTMB(y/N ~ 1, weights=N, data=d, family=binomial())
    expect_equal(fixef(m0),fixef(m3))
    m1 <- glmmTMB((y>5)~1,data=d,family=binomial)
    m2 <- glmmTMB(factor(y>5)~1,data=d,family=binomial)
    expect_equal(c(unname(logLik(m1))),-6.931472,tol=1e-6)
    expect_equal(c(unname(logLik(m2))),-6.931472,tol=1e-6)          

})

test_that("formula expansion", {
              ## test that formulas are expanded in the call/printed
              form <- Reaction ~ Days + (1|Subject)
              expect_equal(grep("Reaction ~ Days",
                       capture.output(print(glmmTMB(form, sleepstudy))),
            fixed=TRUE),1)
})

test_that("NA handling", {
    data(sleepstudy,package="lme4")
    ss <- sleepstudy
    ss$Days[c(2,20,30)] <- NA
    op <- options(na.action=NULL)
    expect_error(glmmTMB(Reaction~Days,ss),"missing values in object")
    op <- options(na.action=na.fail)
    expect_error(glmmTMB(Reaction~Days,ss),"missing values in object")
    expect_equal(unname(fixef(glmmTMB(Reaction~Days,ss,na.action=na.omit))[[1]]),
                 c(249.70505,11.11263),
                 tolerance=1e-6)
    op <- options(na.action=na.omit)
    expect_equal(unname(fixef(glmmTMB(Reaction~Days,ss))[[1]]),
                 c(249.70505,11.11263),
                 tolerance=1e-6)
})

test_that("quine NB fit", {
    quine.nb1 <- MASS::glm.nb(Days ~ Sex/(Age + Eth*Lrn), data = quine)
    quine.nb2 <- glmmTMB(Days ~ Sex/(Age + Eth*Lrn), data = quine,
                         family=nbinom2())
    expect_equal(coef(quine.nb1),fixef(quine.nb2)[["cond"]],
                 tolerance=1e-4)
})
## quine.nb3 <- glmmTMB(Days ~ Sex + (1|Age), data = quine,
##                     family=nbinom2())

test_that("contrasts arg", {
    quine.nb1 <- MASS::glm.nb(Days ~ Sex*Age, data = quine,
                              contrasts=list(Sex="contr.sum",Age="contr.sum"))
    quine.nb2 <- glmmTMB(Days ~ Sex*Age, data = quine,
                         family=nbinom2(),
                         contrasts=list(Sex="contr.sum",Age="contr.sum"))
    expect_equal(coef(quine.nb1),fixef(quine.nb2)[["cond"]],
                 tolerance=1e-4)
})



