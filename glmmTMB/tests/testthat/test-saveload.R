stopifnot(require("testthat"),
          require("glmmTMB"))

context("Saving and loading glmmTMB objects")

test_that("summary consistency", {
    data(sleepstudy, package="lme4")
    fm1 <- glmmTMB(Reaction ~ Days + (1|Subject), sleepstudy)
    s1 <- capture.output(print(summary(fm1)))
    save(fm1, file="fm1.Rdata")
    load("fm1.Rdata")
    file.remove("fm1.Rdata")
    s2 <- capture.output(print(summary(fm1)))
    expect_identical(s1, s2)
})
