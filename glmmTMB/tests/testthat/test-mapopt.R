stopifnot(require("testthat"),
          require("glmmTMB"))

data(Salamanders, package = "glmmTMB")

if (identical(Sys.getenv("NOT_CRAN"), "true")) {
m1 <- glmmTMB(count~ mined, family=poisson, data=Salamanders,
              start=list(beta=c(0,2)),
              map=list(beta=factor(c(1,NA))))

m2 <- glmmTMB(count~ mined + (1|site), family=poisson, data=Salamanders,
              start=list(theta=log(2)),
              map=list(theta=factor(NA)))

m3 <- glmmTMB(count~ mined + (1|site),
              zi = ~1, family=poisson, data=Salamanders,
              start=list(theta=log(2), betazi=c(-1)),
              map=list(theta=factor(NA),
                       betazi=factor(NA)))

m4_nomap <- glmmTMB(count~ mined + (1|site),
              zi=~mined,  family=poisson, data=Salamanders)

m4 <- glmmTMB(count~ mined + (1|site),
              zi=~mined,  family=poisson, data=Salamanders,
              map=list(theta=factor(NA)),
              start = list(theta=log(10)))

data(schizophrenia, package="lme4")
schizophrenia$Week_f <- as.factor(schizophrenia$Week)
schizophrenia$TxDrug_f <- as.factor(schizophrenia$TxDrug)
m5 <- glmmTMB(imps79 ~ Week_f + TxDrug_f:Week_f + (1 | id), data=schizophrenia,
              map=list(beta=factor(c(1:7, NA, 8:13))),
              start=list(beta=rep(0, 14)), # Ensure test works if default args change
              family=Gamma(link = "log"))

m1optim <- update(m1, control=glmmTMBControl(optimizer=optim,
                                             optArgs=list(method="BFGS")))

test_that("basic mapping works", {
    expect_equal(fixef(m1)$cond[[2]], 2.0)
    expect_equal(exp(getME(m2,"theta")), 2.0)
    expect_equal(fixef(m3)$zi[[1]], -1.0)
    expect_equal(fixef(m5)$cond[[8]], 0)
})


test_that("predict works with mapped params", {
          expect_equal(lapply(predict(m1,se.fit=TRUE),unique),
                       list(fit = c(-1.18646939995962, 0.81353060004038), se.fit = 0.0342594326326739),
                       tolerance=1e-6)
          expect_equal(lapply(predict(m5,se.fit=TRUE,newdata=data.frame(Week_f = "0", TxDrug_f = c("0", "1"), id=4720)), unique),
                       list(fit = 1.902812, se.fit = 0.1132216),
                       tolerance =1e-6)
        })

m1_sd <- c(`(Intercept)` = 0.0342594326326741, minedno = NA_real_)

test_that("vcov works with mapped params", {
    expect_equal(dim(vcov(m1, include_nonest=FALSE)$cond),c(1,1))
    expect_equal(dim(vcov(m1, full=TRUE, include_nonest = FALSE)), c(1,1))
    expect_equal(dim(vcov(m2)$cond),c(2,2))
    expect_equal(dim(vcov(m2,full=TRUE, include_nonest=FALSE)),c(2,2))
    v1 <- vcov(m1,include_nonest=TRUE)
    expect_equal(dim(v1$cond),c(2,2))
    expect_equal(sqrt(diag(v1$cond)), m1_sd, tolerance=1e-6)
    v5 <- vcov(m5,include_nonest=TRUE)
    expect_equal(dim(v5$cond),c(14, 14))
    expect_equal(dim(vcov(m5, include_nonest=FALSE)$cond),c(13, 13))
    expect_all_equal(v5$cond[8, ],NA_real_)
    expect_all_equal(v5$cond[, 8],NA_real_)
})

test_that("summary works with mapped params", {
    expect_equal(summary(m1)$coef$cond[,"Std. Error"], m1_sd)
    expect_equal(coef(summary(m5))$cond[8,],
                 c(Estimate = 0, `Std. Error` = NA, `z value` = NA, `Pr(>|z|)` = NA))
})

test_that("confint works with mapped params", {
    cm1 <- confint(m1)
    expect_equal(dim(cm1), c(1,3))
    expect_equal(rownames(cm1), "(Intercept)")
    cm2 <- confint(m2)
    expect_equal(dim(cm2), c(2,3))
    expect_equal(rownames(cm2), c("(Intercept)","minedno"))
    cm3 <- confint(m3)
    expect_equal(dim(cm3), c(2,3))
    expect_equal(rownames(cm3), c("(Intercept)","minedno"))
    cm4 <- confint(m4)
    expect_equal(dim(cm4), c(4,3))
    expect_equal(rownames(cm4),
       c("cond.(Intercept)", "cond.minedno", "zi.(Intercept)", "zi.minedno"))
    cm4_nomap <- confint(m4_nomap)
    cm5 <- confint(m5, include_nonest = TRUE)
    expect_equal(dim(cm5), c(15,3))
    expect_equal(cm5[8, ], c(`2.5 %` = NA, `97.5 %` = NA, Estimate = 0))
})


context("alternate optimizers")

test_that("alternate optimizers work", {
          expect_equal(fixef(m1),fixef(m1optim),
                       tol=1e-4)
          expect_false(identical(fixef(m1),fixef(m1optim)))
})

test_that("summary", {
    expect_equal(coef(summary(m1))$cond["minedno",],
                 c(Estimate = 2, `Std. Error` = NA, `z value` = NA, `Pr(>|z|)` = NA))
    expect_equal(coef(summary(m3))$zi["(Intercept)",],
                 c(Estimate = -1, `Std. Error` = NA, `z value` = NA, `Pr(>|z|)` = NA))

})

} ## skip on CRAN
