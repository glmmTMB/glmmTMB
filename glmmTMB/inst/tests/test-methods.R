stopifnot(require("testthat"),
          require("glmmTMB"))

data(sleepstudy, cbpp,
     package = "lme4")

context("Basic methods tests")

fm2 <- glmmTMB(Reaction ~ Days + (Days| Subject), sleepstudy)
fm2P <- glmmTMB(round(Reaction) ~ Days + (Days| Subject), sleepstudy,
               family=poisson)
fm2NB <- glmmTMB(round(Reaction) ~ Days + (Days| Subject), sleepstudy,
               family=list(family="nbinom2",link="log"))
               
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
          })


test_that("VarCorr", {
   vv <- VarCorr(fm2)
   vv2 <- vv$cond$Subject
   expect_equal(dim(vv2),c(2,2))
   expect_equal(outer(attr(vv2,"stddev"),attr(vv2,"stddev"))*attr(vv2,"correlation"),
                vv2,check.attributes=FALSE)
}
