## check that everything works in weird environments

stopifnot(require("testthat"),
          require("glmmTMB"))

data(sleepstudy, cbpp,
     package = "lme4")

## need global env for test_that
sleepstudy <<- transform(sleepstudy, DaysFac = factor(Days))

context("basic examples")

test_that("basic example #1", {
    fitFun <- function(dat){
        glmmTMB(Reaction ~ Days + (1|Subject), data=dat)
    }
    f0 <- glmmTMB(Reaction ~ Days + (1|Subject), data=sleepstudy)
    f1 <- fitFun(sleepstudy)
    uncall <- function(x) {
        x$call <- NULL
        return(x)
    }
    expect_equal(uncall(f0),uncall(f1))
})

test_that("paranoia", {
    formFun <- function() {
        return(Reaction ~ Days + (1|Subject))
    }
    fitFun <- function(f,dat){
        glmmTMB(f, data=dat)
    }
    f0 <- glmmTMB(Reaction ~ Days + (1|Subject), data=sleepstudy)
    f1 <- fitFun(formFun(),sleepstudy)
    uncall <- function(x) {
        x$call <- NULL
        return(x)
    }
    expect_equal(uncall(f0),uncall(f1))
})
test_that("dispformula env", {
	fitFun2 <- function(dat){
        glmmTMB(count ~ 1, data=dat, family="poisson" )
    }
    m0 <- fitFun2(Salamanders)
    m1 <- glmmTMB(count ~ 1, data= Salamanders, family="poisson")
    uncall <- function(x) {
        x$call <- NULL
        return(x)
    }
    expect_equal(uncall(summary(m0)), uncall(summary(m1)))
})
