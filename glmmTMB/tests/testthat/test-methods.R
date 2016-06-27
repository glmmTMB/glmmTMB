stopifnot(require("testthat"),
          require("glmmTMB"))

data(sleepstudy, cbpp,
     package = "lme4")

fm2 <- glmmTMB(Reaction ~ Days + (Days| Subject), sleepstudy)
fm0 <- update(fm2, . ~ . -Days)
fm2P <- glmmTMB(round(Reaction) ~ Days + (Days| Subject), sleepstudy,
               family=poisson)

context("basic methods")

## gives warnings (crazy model ...)
fm2NB <- suppressWarnings(
    glmmTMB(round(Reaction) ~ Days + (Days| Subject), sleepstudy,
               family=list(family="nbinom2",link="log")))
               
test_that("Fitted and residuals", {
    expect_equal(length(fitted(fm2)),nrow(sleepstudy))
    expect_equal(mean(fitted(fm2)),298.507891)
    expect_equal(mean(residuals(fm2)),0,tol=1e-5)
    ## Pearson and response are the same for a Gaussian model
    expect_equal(residuals(fm2,type="response"),
                 residuals(fm2,type="pearson"))
    ## ... but not for Poisson ...
    expect_false(mean(residuals(fm2P,type="response"))==
                 mean(residuals(fm2P,type="pearson")))
    expect_error(residuals(fm2NB,type="pearson"),
                 "variance function undefined")
})

test_that("Predict", {
     expect_equal(predict(fm2),predict(fm2,newdata=sleepstudy))
 })


test_that("VarCorr", {
   vv <- VarCorr(fm2)
   vv2 <- vv$cond$Subject
   expect_equal(dim(vv2),c(2,2))
   expect_equal(outer(attr(vv2,"stddev"),
                      attr(vv2,"stddev"))*attr(vv2,"correlation"),
                vv2,check.attributes=FALSE)
})

test_that("drop1", {
      dd <- drop1(fm2,test="Chisq")
      expect_equal(dd$AIC,c(1763.94,1785.48),tol=1e-4)              
          })
test_that("anova", {
      aa <- anova(fm0,fm2)
      expect_equal(aa$AIC,c(1785.48,1763.94),tol=1e-4)
          })

test_that("terms", {
    ## test whether terms() are returned with predvars for doing
    ## model prediction etc. with complex bases
    dd <<- data.frame(x=1:10,y=1:10)
    require("splines")
    m <- glmmTMB(y~ns(x,3),dd)
    ## if predvars is not properly attached to term, this will
    ## fail as it tries to construct a 3-knot spline from a single point
    model.matrix(delete.response(terms(m)),data=data.frame(x=1))
})

