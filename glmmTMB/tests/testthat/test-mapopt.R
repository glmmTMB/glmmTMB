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

m1optim <- update(m1, control=glmmTMBControl(optimizer=optim,
                                             optArgs=list(method="BFGS")))

test_that("basic mapping works", {
    expect_equal(fixef(m1)$cond[[2]], 2.0)
    expect_equal(exp(getME(m2,"theta")), 2.0)
})

test_that("predict works with mapped params",
          expect_equal(vapply(predict(m1,se.fit=TRUE),unique,numeric(1)),
                       c(fit = -1.18646939995962, se.fit = 0.0342594326326737),
                       tolerance=1e-6))

context("alternate optimizers")

test_that("alternate optimizers work", {
          expect_equal(fixef(m1),fixef(m1optim),
                       tol=1e-4)
          expect_false(identical(fixef(m1),fixef(m1optim)))
})
          

        
