stopifnot(require("testthat"),
          require("glmmTMB"))

data(Salamanders, package = "glmmTMB")
context("mapping")

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

m1optim <- update(m1, control=glmmTMBControl(optimizer=optim,
                                             optArgs=list(method="BFGS")))

test_that("basic mapping works", {
    expect_equal(fixef(m1)$cond[[2]], 2.0)
    expect_equal(exp(getME(m2,"theta")), 2.0)
    expect_equal(fixef(m3)$zi[[1]], -1.0)
})

test_that("predict works with mapped params",
          expect_equal(vapply(predict(m1,se.fit=TRUE),unique,numeric(1)),
                       c(fit = -1.18646939995962, se.fit = 0.0342594326326737),
                       tolerance=1e-6)
          )

test_that("vcov works with mapped params", {
    expect_equal(dim(vcov(m1)$cond),c(1,1))
    expect_equal(dim(vcov(m1,full=TRUE)),c(1,1))
    expect_equal(dim(vcov(m2)$cond),c(2,2))
    expect_equal(dim(vcov(m2,full=TRUE)),c(2,2))
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
})


context("alternate optimizers")

test_that("alternate optimizers work", {
          expect_equal(fixef(m1),fixef(m1optim),
                       tol=1e-4)
          expect_false(identical(fixef(m1),fixef(m1optim)))
})


        
