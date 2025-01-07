stopifnot(require("testthat"), require("glmmTMB"))

m1 <- glmmTMB(mpg ~ 1 + hp + disp + drat + wt + qsec,
        family = "gaussian",
        data = mtcars,
        map = list(beta = factor(c(1, 2, 2, 3, 3, 3))))

test_that("mapping equal coefs works", {
    expect_identical(fixef(m1)$cond[[2]], fixef(m1)$cond[[3]])
    expect_identical(fixef(m1)$cond[[4]], fixef(m1)$cond[[5]])
    expect_identical(fixef(m1)$cond[[5]], fixef(m1)$cond[[6]])
})

v <- vcov(m1)
cc <- cov2cor(v$cond)
test_that("vcov for equal-mapping", {
    expect_equal(cc[2,3], 1.0)
    expect_true(all(abs(cc[4:6, 4:6] - 1) < 1e-8))
    expect_equal(v$cond[2,2], v$cond[3,3])
    expect_equal(v$cond[4,4], v$cond[5,5])
    expect_equal(v$cond[5,5], v$cond[6,6])
})

vcov(m1, full = TRUE)
