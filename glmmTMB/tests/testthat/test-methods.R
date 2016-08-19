stopifnot(require("testthat"),
          require("glmmTMB"))

data(sleepstudy, cbpp,
     package = "lme4")

## FIXME: fit these centrally and restore, to save time
fm2   <- glmmTMB(Reaction ~ Days + (Days| Subject), sleepstudy)
fm0   <- update(fm2, . ~ . -Days)
fm2P  <- update(fm2, round(Reaction) ~ ., family=poisson)
fm2G  <- update(fm2, family=Gamma(link="log"))
fm2NB <- update(fm2P, family=nbinom2)
## for testing sigma() against base R
fm3   <- update(fm2, . ~ Days)
fm3G <-  update(fm3, family=Gamma(link="log"))
fm3NB <- update(fm3, round(Reaction) ~ ., family=nbinom2)

context("basic methods")

test_that("Fitted and residuals", {
    expect_equal(length(fitted(fm2)),nrow(sleepstudy))
    expect_equal(mean(fitted(fm2)),298.507891)
    expect_equal(mean(residuals(fm2)),0,tol=1e-5)
    ## Pearson and response are the same for a Gaussian model
    expect_equal(residuals(fm2,type="response"),
                 residuals(fm2,type="pearson"))
    ## ... but not for Poisson or NB ...
    expect_false(mean(residuals(fm2P,type="response"))==
                 mean(residuals(fm2P,type="pearson")))
    expect_false(mean(residuals(fm2NB,type="response"))==
                 mean(residuals(fm2NB,type="pearson")))
})

test_that("Predict", {
    expect_equal(predict(fm2),predict(fm2,newdata=sleepstudy))
    pr2se <- predict(fm2, se.fit=TRUE)
    i <- sample(nrow(sleepstudy), 20)
    newdata <- sleepstudy[i, ]
    pr2sub <- predict(fm2, newdata=newdata, se.fit=TRUE)
    expect_equivalent(pr2se$fit, predict(fm2))
    expect_equivalent(pr2se$fit[i], pr2sub$fit)
    expect_equivalent(pr2se$se.fit[i], pr2sub$se.fit)
    expect_equal(unname( pr2se$   fit[1] ), 254.2208, tol=1e-4)
    expect_equal(unname( pr2se$se.fit[1] ), 12.94514, tol=1e-4)
    expect_equal(unname( pr2se$   fit[100] ), 457.9684, tol=1e-4)
    expect_equal(unname( pr2se$se.fit[100] ), 14.13943, tol=1e-4)
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

test_that("summary_print", {
    ## no dispersion printed for Gaussian or disp==1 families
    c1 <<- capture.output(print(summary(fm2)))
    expect_true(!any(grepl("Dispersion",c1)))
    c2 <<- capture.output(print(summary(fm2P)))
    expect_true(!any(grepl("Dispersion",c2)))
    ## extract numeric dispersion from printed output 
    c3 <<- capture.output(print(summary(fm2G)))
    dline <- grep("Dispersion",c3,value=TRUE)
    expect_equal(as.numeric(gsub("[^0-9.]","",dline)),0.00853)
})

test_that("sigma", {
    s1 <<- sigma(lm(Reaction~Days,sleepstudy))
    s2 <<- sigma(glm(Reaction~Days,sleepstudy,family=Gamma(link="log")))
    s3 <<- MASS::glm.nb(round(Reaction)~Days,sleepstudy)
    ## remove bias-correction
    expect_equal(sigma(fm3),s1*(1-1/nobs(fm3)),tolerance=1e-3)
    expect_equal(sigma(fm3G),s2,tolerance=5e-3)
    expect_equal(s3$theta,sigma(fm3NB),tolerance=1e-4)
})
