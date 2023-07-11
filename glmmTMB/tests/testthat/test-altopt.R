stopifnot(require("testthat"),
          require("glmmTMB"))

data(sleepstudy, cbpp,
     package = "lme4")


context("alternative optimizers")

test_that("downstream methods work with optim()", {
  skip_on_cran()
  m1 <- glmmTMB(count~ mined, family = poisson, data = Salamanders,
                control = glmmTMBControl(optimizer = optim,
                                         optArgs = list(method="BFGS")))
    expect_is(summary(m1),"summary.glmmTMB")
    expect_is(confint(m1),"matrix")
})
