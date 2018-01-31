stopifnot(require("testthat"),
          require("glmmTMB"))

context("offsets")

set.seed(101)
n <- 10000
mux <- 10
sdx <- 10
a <- .1
b <- 0
residsd <- .01
x <- rnorm(n, mux, sdx)
o <<- 2*x+100
o2 <- rep(c(1,5),each=n/2)
y <- a*x+b+o+rnorm(n, sd=residsd)
y2 <- a*x+b+rnorm(n, sd=residsd*o2)
## global assignment for testthat
dat <<- data.frame(y, y2, x, o, o2)
m.lm <- lm(y~x, offset=o, dat)

test_that("LM with offset as argument", {
    m1 <- glmmTMB(y~x, offset=o, dat) 
    expect_equal(fixef(m1)[[1]], coef(m.lm), tol=1e-6)
    m3 <- glmmTMB(y~x, offset=o)
    expect_equal(fixef(m3)[[1]], coef(m.lm), tol=1e-6)
})

test_that("LM with offset in formula", {
    m2 <- glmmTMB(y~x+offset(o), dat)
    expect_equal(fixef(m2)[[1]], coef(m.lm), tol=1e-6)
    m4 <- glmmTMB(y~x+offset(o))
    expect_equal(fixef(m4)[[1]], coef(m.lm), tol=1e-6)
})

test_that("LM with offset in zero-inflation formula", {
    ## don't have anything sensible to try here yet ...
    ## glmmTMB(y~x,zi=~1+offset(o), dat)
})

test_that("LM with offset in dispersion formula", {
    expect_equal(sigma(glmmTMB(y~x, dat)),
                 sigma(glmmTMB(y2~x,disp=~1+offset(log(o2)*2), dat)),
                 tolerance=1e-3)
})
