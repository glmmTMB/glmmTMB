stopifnot(require("testthat"), require("glmmTMB"), require("mgcv"))

smooth_start_data <- function() {
    set.seed(928)
    d <- data.frame(x = runif(240), z = runif(240),
                    g = factor(rep(seq_len(12), 20)))
    d$y <- 100 + 20*sin(6*d$x) + 10*cos(6*d$z) +
        rnorm(12, sd=5)[d$g] + rnorm(nrow(d), sd=3)
    d
}

test_that("only smooth covariance starts change, in formula order", {
    d <- smooth_start_data()
    for (f in list(y ~ (1+x|g) + s(x,k=6) + s(z,k=5),
                   y ~ s(x,k=6) + (1+x|g) + s(z,k=5),
                   y ~ s(x,k=6) + s(z,k=5) + (1+x|g))) {
        a <- glmmTMB(f, data=d, doFit=FALSE)
        info <- a$condList$reTrms$smooth_info
        positions <- which(lengths(info) > 0L)
        sizes <- vapply(a$condReStruc, `[[`, numeric(1), "blockNumTheta")
        smooth_theta <- cumsum(c(0, sizes))[positions] + 1L
        expect_true(all(is.finite(a$parameters$theta[smooth_theta])))
        expect_true(all(a$parameters$theta[smooth_theta] != 0))
        expect_equal(a$parameters$theta[-smooth_theta], rep(0, 3))
        manual <- glmmTMB(f, data=d, doFit=FALSE,
                         start=list(theta=rep(0, sum(sizes))))
        expect_identical(a$parameters[names(a$parameters) != "theta"],
                         manual$parameters[names(manual$parameters) != "theta"])
        expect_identical(manual$parameters$theta, rep(0, sum(sizes)))
        expect_identical(a$parameters$betadisp, 0)
    }
})

test_that("explicit starts and mapped parameters retain their meaning", {
    d <- smooth_start_data()
    f <- y ~ s(x,k=6) + (1|g)
    base <- glmmTMB(f, data=d, doFit=FALSE)
    supplied <- list(beta=c(7, 8), b=rep(.4, ncol(base$data.tmb$Z)),
                     theta=c(-2, .3), betadisp=2)
    a <- glmmTMB(f, data=d, doFit=FALSE, start=supplied)
    expect_identical(a$parameters[names(supplied)], supplied)
    a <- glmmTMB(f, data=d, doFit=FALSE, map=list(theta=factor(c(NA, 1))))
    expect_identical(a$parameters$theta, c(0, 0))
})

test_that("initialization follows response units and fitting rows", {
    d <- smooth_start_data()
    f <- y ~ s(x,k=6)
    a <- glmmTMB(f, data=d, doFit=FALSE)
    scaled <- transform(d, y=10*y + 300)
    b <- glmmTMB(f, data=scaled, doFit=FALSE)
    expect_equal(b$parameters$theta - a$parameters$theta, log(10))
    expect_identical(b$parameters$betadisp, a$parameters$betadisp)

    d$w <- seq(.5, 1.5, length.out=nrow(d))
    d$x[3] <- NA
    subset <- d$z < .8
    a <- glmmTMB(f, data=d, weights=w, subset=subset, doFit=FALSE,
                 na.action=na.omit)
    used <- d[subset & !is.na(d$x), ]
    b <- glmmTMB(f, data=used, weights=w, doFit=FALSE)
    expect_equal(a$parameters, b$parameters)
    d <- used
    d$off <- seq(-2, 2, length.out=nrow(d))
    a <- glmmTMB(y ~ s(x,k=6) + offset(off), data=d, doFit=FALSE)
    b <- glmmTMB(f, data=transform(d, y=y-off), doFit=FALSE)
    expect_equal(a$parameters, b$parameters)
})

test_that("models outside the prototype retain their default starts", {
    d <- smooth_start_data()
    d$y <- round(d$y)
    cases <- list(
        list(y ~ x + (1|g), family=gaussian()),
        list(y ~ s(x,k=6), family=gaussian(link="log")),
        list(y ~ s(x,k=6), family=poisson(link="sqrt")),
        list(y ~ s(x,k=6), family=Gamma(link="identity")),
        list(y ~ s(x,k=6), family=genpois()),
        list(y ~ s(x,k=6), dispformula=~z),
        list(y ~ s(x,k=6), dispformula=~0),
        list(y ~ s(x,k=6), ziformula=~1))
    for (case in cases) {
        a <- do.call(glmmTMB, c(case, list(data=d, doFit=FALSE)))
        expect_true(all(a$parameters$theta == 0))
        expect_identical(glmmTMB:::smoothStart(a$parameters, a$condList,
            a$condReStruc, a$data.tmb, a$family, a$mapArg), a$parameters)
    }
    d$y <- 5
    a <- glmmTMB(y ~ s(x,k=6), data=d, doFit=FALSE)
    expect_identical(a$parameters$theta, 0)
})

test_that("volcano smooths fit without explicit starts under ML and REML", {
    d <- data.frame(z=as.vector(volcano), x=as.vector(row(volcano)),
                    y=as.vector(col(volcano)))
    set.seed(1)
    d$z <- d$z + rnorm(length(volcano), sd=15)
    d <- d[sample(nrow(d), 100), ]
    for (reml in c(FALSE, TRUE)) {
        a <- glmmTMB(z ~ s(x,y), data=d, REML=reml)
        b <- mgcv::gam(z ~ s(x,y), data=d, method=if (reml) "REML" else "ML")
        expect_equal(a$fit$convergence, 0)
        expect_true(a$sdr$pdHess)
        expect_lt(max(abs(predict(a) - predict(b))), 1e-3)
        expect_equal(predict(a, newdata=d), predict(a), tolerance=1e-7)
    }
    explicit <- glmmTMB(z ~ s(x,y,bs="tp"), data=d)
    expect_true(explicit$sdr$pdHess)
})

test_that("transformed sleepstudy smooth fits without explicit starts", {
    data("sleepstudy", package="lme4")
    a <- glmmTMB(Reaction ~ s(log(Days+1),k=4), data=sleepstudy, REML=TRUE)
    b <- mgcv::gam(Reaction ~ s(log(Days+1),k=4), data=sleepstudy, method="REML")
    expect_equal(a$fit$convergence, 0)
    expect_true(a$sdr$pdHess)
    expect_lt(max(abs(predict(a) - predict(b))), 1e-3)
    expect_equal(-as.numeric(logLik(a)), as.numeric(b$gcv.ubre), tolerance=1e-7)
    expect_equal(predict(a, newdata=sleepstudy), predict(a), tolerance=1e-7)
})
