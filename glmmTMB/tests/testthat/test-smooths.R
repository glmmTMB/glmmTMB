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
    expect_silent(g1 <- glmmTMB(x ~ 1 + s(y, x, bs = "sos"),
                  dispformula = ~0, data = dd))

})

test_that("smooth + diag() specials", {
    set.seed(101)
    n <- 100
    dd <- data.frame(y = runif(n, -80, 80), x = runif(n, -180, 180),
                     f1 = factor(sample(1:4, size = n, replace = TRUE)),
                     f2 = factor(sample(1:4, size = n, replace = TRUE)))
    ff <- z ~ 1 + s(y, x, bs = "sos") + diag(f1|f2)
    dd$z <- simulate_new(ff[-2],
                         dispformula = ~ 0,
                         newdata = dd,
                         newparams = list(beta = 1, theta = rep(1, 5)),
                         family = gaussian)[[1]]

    g1 <- glmmTMB(ff, dispformula = ~0, data = dd)
    expect_equal(getME(g1, "theta"),
                 c(0.829644130428793, 0.789420770119563,
                   1.1073848292829, 0.522906430406326, 
                   1.16293701025298),
                 tolerance = 1e-4)
})

test_that("multiple smooths", {

    fake_data <- data.frame(y = rnorm(n = 100),
                            x1 = rnorm(n = 100),
                            x2 = rnorm(n = 100))
    
    g1 <- glmmTMB(
        y ~ s(x1) + s(x2),
        data = fake_data
    )
    expect_equal(names(fixef(g1)$cond),
                       c("(Intercept)", "s(x1)1", "s(x2)1"))

})
