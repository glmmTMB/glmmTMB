stopifnot(require("testthat"),
          require("glmmTMB"),
          require("mgcv"))

data("Nile")
ndat <- data.frame(time = time(Nile), val = c(Nile))

test_that("basic smooth", {
    sm0 <- gam(val ~ s(time), data = ndat, method = "REML")
    ## debug(glmmTMB:::getXReTrms)
    sm1 <- glmmTMB(val ~ s(time), data = ndat,
                   REML = TRUE, start = list(theta = log(sd(ndat$val))))
    expect_equal(predict(sm1), predict(sm1, newdata = ndat))
    expect_equal(unname(c(predict(sm0))), predict(sm1), tolerance = 5e-5)
    expect_equal(-1*c(logLik(sm1)), unname(c(sm0$gcv.ubre)))
    r1 <- ranef(sm1)$cond[[1]]
    r2 <- split(coef(sm0), substr(names(coef(sm0)), 1, 3))[[2]]
    ## not finished ...
    ## expect_equal(unname(unlist(r1[[1]])), unname(r2[[1]]), tolerance = 1e-5)
})

test_that("smooth with no fixed-effect components", {
    set.seed(101)
    n <- 100
    dd <- data.frame(y = runif(n, -80, 80), x = runif(n, -180, 180))
    dd$z <- simulate_new(~ 1 + s(y, x, bs = "sos"),
                         dispformula = ~ 0,
                         newdata = dd,
                         newparams = list(beta = 1, theta = 1),
                         family = gaussian)[[1]]
    ## fails if n==20
    ## Error in cbind(S, matrix(0, k - length(ind), length(ind))) : 
    ## number of rows of matrices must match (see arg 2)

    ## fits without warnings but not sensible ...
    g1 <- glmmTMB(x ~ 1 + s(y, x, bs = "sos"),
                  dispformula = ~0, data = dd)

})
