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
r1 <- rnorm(n, sd=residsd)
r2 <- rnorm(n, sd=residsd*o2)
y0 <- a*x + b + r1
y1 <- a*x+b+o+r1
y2 <- a*x+b+r2
y3 <- a*x + b + o + r2
## global assignment for testthat
dat <<- data.frame(y0, y1, y2, y3, x, o, o2, o3=o)
m.lm <- lm(y1~x, offset=o, dat)
m.lm0 <- lm(y1~x, dat)

test_that("LM with offset as argument", {
    m1 <- glmmTMB(y1~x, offset=o, dat) 
    expect_equal(fixef(m1)[[1]], coef(m.lm), tol=1e-6)
    m3 <- glmmTMB(y1~x, offset=o)
    expect_equal(fixef(m3)[[1]], coef(m.lm), tol=1e-6)
})

test_that("LM with offset in formula", {
    m2 <- glmmTMB(y1~x+offset(o), dat)
    expect_equal(fixef(m2)[[1]], coef(m.lm), tol=1e-6)
    m4 <- glmmTMB(y1~x+offset(o))
    expect_equal(fixef(m4)[[1]], coef(m.lm), tol=1e-6)
})

## test_that("LM with offset in zero-inflation formula", {
    ## don't have anything sensible to try here yet ...
    ## glmmTMB(y~x,zi=~1+offset(o), dat)
## })

test_that("LM with offset in formula - variable not in environment", {
    m5 <- glmmTMB(y1~x,offset=o3, dat)
    expect_equal(fixef(m5)[[1]],coef(m.lm), tol=1e-6)
    nullvalue <- NULL
    m6 <- glmmTMB(y1~x,offset=nullvalue, dat)
    expect_equal(fixef(m6)[[1]],coef(m.lm0), tol=1e-6)
})

test_that("LM with offset in dispersion formula", {
    expect_equal(sigma(glmmTMB(y1~x, dat)),
                 sigma(glmmTMB(y2~x,disp=~1+offset(log(o2)*2), dat)),
                 tolerance=1e-3)
    
})

test_that("LM with multiple offsets (cond/dispersion)", {
    m1 <<- glmmTMB(y0~x, dat)
    m2 <<- glmmTMB(y3~x+offset(o),disp=~1+offset(log(o2)*2), dat)
    expect_equal(sigma(m1),sigma(m2),tolerance=1e-3)
    expect_equal(fixef(m1),fixef(m2),tolerance=1e-3)
})

## this was broken by an earlier multiple-offset formulation
test_that("LM with random crap in the formula", {
    m1 <<- glmmTMB(y0~dat$x, dat)
    m2 <<- glmmTMB(y0~x, dat)
    expect_equal(unname(fixef(m1)$cond),unname(fixef(m2)$cond))
})

