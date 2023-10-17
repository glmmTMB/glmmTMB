stopifnot(require("testthat"),
          require("glmmTMB"),
          require("mgcv"))

data("Nile")
ndat <- data.frame(time = time(Nile), val = c(Nile))

test_that("basic smooth", {
    sm0 <- gam(val ~ s(time), data = ndat, method = "REML")
    sm1 <- glmmTMB(val ~ s(time), data = ndat,
                   REML = TRUE, start = list(theta = 5))
    expect_equal(predict(sm1), predict(sm1, newdata = ndat))
    expect_equal(unname(c(predict(sm0))), predict(sm1))
    expect_equal(-1*c(logLik(sm1)), unname(c(sm0$gcv.ubre)))
    r1 <- ranef(sm1)$cond[[1]]
    r2 <- split(coef(sm0), substr(names(coef(sm0)), 1, 3))[[2]]
    ## not finished ...
    ## expect_equal(unname(unlist(r1[[1]])), unname(r2[[1]]), tolerance = 1e-5)
})
