stopifnot(require("testthat"), require("glmmTMB"), require("mgcv"))

smooth_start_glm_data <- function(case, dimension=1L, seed=1L) {
    set.seed(928+seed)
    d <- data.frame(x=runif(400), z=runif(400))
    signal <- if (dimension == 1L) 1.5*sin(2*pi*d$x) else
        1.5*sin(2*pi*d$x)*cos(2*pi*d$z)
    eta <- signal + switch(case, poisson_large=log(1000),
                           bernoulli_rare=-3.5, bernoulli=0,
                           binomial=0, beta=0, 1)
    mu <- if (case %in% c("bernoulli", "bernoulli_rare", "binomial", "beta"))
        plogis(eta) else exp(eta)
    family <- switch(case, poisson=, poisson_large=poisson(),
                     bernoulli=, bernoulli_rare=, binomial=binomial(),
                     nbinom1=nbinom1(), nbinom2=nbinom2(),
                     Gamma_log=Gamma(link="log"), Gamma_inverse=Gamma(),
                     beta=beta_family(), tweedie=tweedie())
    trials <- if (case == "binomial") 20 else 1
    d$y <- switch(case,
        poisson=, poisson_large=rpois(nrow(d),mu),
        bernoulli=, bernoulli_rare=, binomial=rbinom(nrow(d),trials,mu),
        nbinom1=rnbinom(nrow(d),mu=mu,size=mu/2),
        nbinom2=rnbinom(nrow(d),mu=mu,size=2),
        Gamma_log=, Gamma_inverse=rgamma(nrow(d),shape=3,scale=mu/3),
        beta=rbeta(nrow(d),mu*20,(1-mu)*20),
        tweedie=mgcv::rTweedie(mu,p=1.5,phi=.7))
    if (case == "binomial") d$failures <- trials-d$y
    formula <- if (dimension == 1L) y ~ s(x,k=8) else y ~ s(x,z,k=20)
    if (case == "binomial") formula[[2L]] <- quote(cbind(y,failures))
    list(formula=formula,data=d,family=family)
}

test_that("GLM initialization changes only conditional smooth covariance starts", {
    for (case in c("poisson", "bernoulli", "binomial", "nbinom1", "nbinom2",
                   "Gamma_log", "Gamma_inverse", "beta", "tweedie")) {
        args <- smooth_start_glm_data(case)
        a <- do.call(glmmTMB,c(args,list(doFit=FALSE)))
        old <- do.call(glmmTMB,c(args,list(doFit=FALSE,start=list(theta=0))))
        expect_true(all(is.finite(a$parameters$theta)))
        expect_false(identical(a$parameters$theta,0))
        expect_identical(a$parameters[names(a$parameters)!="theta"],
                         old$parameters[names(old$parameters)!="theta"])
        expect_identical(old$parameters$theta,0)
        mapped <- do.call(glmmTMB,c(args,list(doFit=FALSE,
                                   map=list(theta=factor(NA)))))
        expect_identical(mapped$parameters,old$parameters)
    }
})

test_that("GLM fits match the former defaults with 1D and 2D smooths", {
    skip_on_cran()
    cases <- c("poisson", "poisson_large", "bernoulli", "bernoulli_rare",
               "binomial", "nbinom1", "nbinom2", "Gamma_log", "beta", "tweedie")
    for (case in cases) for (dimension in 1:2) for (seed in 1:3) {
        args <- smooth_start_glm_data(case,dimension,seed)
        a <- do.call(glmmTMB,args)
        old <- do.call(glmmTMB,c(args,list(start=list(theta=0))))
        info <- paste(case,dimension,seed)
        expect_equal(a$fit$convergence,0,info=info)
        expect_true(a$sdr$pdHess,info=info)
        expect_true(abs(a$fit$objective-old$fit$objective)<1e-4,info=info)
    }
})

test_that("binomial starts respect equivalent response encodings", {
    args <- smooth_start_glm_data("binomial")
    d <- args$data
    a <- glmmTMB(cbind(y,failures) ~ s(x,k=8),family=binomial(),
                 data=d,doFit=FALSE)
    d$n <- d$y+d$failures
    d$proportion <- d$y/d$n
    b <- glmmTMB(proportion ~ s(x,k=8),weights=n,family=binomial(),
                 data=d,doFit=FALSE)
    expect_identical(a$parameters,b$parameters)
})

test_that("log-link offsets leave other starts unchanged", {
    args <- smooth_start_glm_data("poisson")
    d <- args$data
    d$exposure <- seq(.5,2,length.out=nrow(d))
    set.seed(929)
    d$y <- rpois(nrow(d),d$exposure*exp(1+1.5*sin(2*pi*d$x)))
    f <- y ~ s(x,k=8)+offset(log(exposure))
    a <- glmmTMB(f,family=poisson(),data=d,doFit=FALSE)
    old <- glmmTMB(f,family=poisson(),data=d,doFit=FALSE,start=list(theta=0))
    expect_identical(a$parameters[names(a$parameters)!="theta"],
                     old$parameters[names(old$parameters)!="theta"])
    ## A link-scale offset must not be subtracted from the response when
    ## computing the response-scale variance heuristic.
    no_offset <- glmmTMB(y ~ s(x,k=8),family=poisson(),data=d,doFit=FALSE)
    expect_identical(a$parameters$theta,no_offset$parameters$theta)
    fit <- glmmTMB(f,family=poisson(),data=d)
    former <- glmmTMB(f,family=poisson(),data=d,start=list(theta=0))
    expect_equal(fit$fit$convergence,0)
    expect_true(fit$sdr$pdHess)
    expect_lt(abs(fit$fit$objective-former$fit$objective),1e-4)
})

test_that("GLM component and ordinary random-effect starts are preserved", {
    args <- smooth_start_glm_data("nbinom2")
    d <- args$data
    d$g <- factor(rep(seq_len(20),20))
    a <- glmmTMB(y ~ (1+x|g) + s(x,k=6),ziformula=~s(z,k=5),
                 dispformula=~s(z,k=5),family=nbinom2(),data=d,doFit=FALSE)
    old <- glmmTMB(y ~ (1+x|g) + s(x,k=6),ziformula=~s(z,k=5),
                   dispformula=~s(z,k=5),family=nbinom2(),data=d,doFit=FALSE,
                   start=list(theta=rep(0,4)))
    expect_identical(a$parameters$theta[1:3],rep(0,3))
    expect_false(identical(a$parameters$theta[4],0))
    expect_identical(a$parameters[names(a$parameters)!="theta"],
                     old$parameters[names(old$parameters)!="theta"])
    supplied <- lapply(old$parameters,function(x) rep(.2,length(x)))
    manual <- glmmTMB(y ~ (1+x|g) + s(x,k=6),ziformula=~s(z,k=5),
                      dispformula=~s(z,k=5),family=nbinom2(),data=d,doFit=FALSE,
                      start=supplied)
    expect_identical(manual$parameters,supplied)
})

test_that("non-Gaussian no-smooth models retain all former starts", {
    for (case in c("poisson", "bernoulli", "binomial", "nbinom1", "nbinom2",
                   "Gamma_log", "Gamma_inverse", "beta", "tweedie")) {
        args <- smooth_start_glm_data(case)
        args$formula[[3L]] <- quote(x)
        a <- do.call(glmmTMB,c(args,list(doFit=FALSE)))
        expect_identical(glmmTMB:::smoothStart(a$parameters,a$condList,
            a$condReStruc,a$data.tmb,a$family,a$mapArg),a$parameters)
    }
})

test_that("GLM fits with zero inflation and varying dispersion retain their optima", {
    skip_on_cran()
    for (case in c("zinb", "nbinom_dispersion", "beta_dispersion", "poisson_slopes")) {
        set.seed(1006)
        d <- data.frame(x=runif(1200),z=runif(1200),
                        g=factor(rep(seq_len(20),60)))
        b <- rnorm(20,sd=.45)
        args <- list(formula=y ~ s(x,k=6)+(1|g))
        if (case == "beta_dispersion") {
            mu <- plogis(.4+sin(6*d$x)+b[d$g])
            phi <- exp(2+.6*cos(6*d$z))
            d$y <- rbeta(nrow(d),mu*phi,(1-mu)*phi)
            args$family <- beta_family()
            args$dispformula <- ~s(z,k=5)
        } else if (case == "poisson_slopes") {
            slopes <- rnorm(20,sd=.4)
            d$y <- rpois(nrow(d),exp(1+.6*sin(6*d$x)+b[d$g]+slopes[d$g]*d$z))
            args$formula <- y ~ s(x,k=6)+(1+z|g)
            args$family <- poisson()
        } else {
            phi <- exp(.7+.6*cos(6*if (case == "zinb") d$x else d$z))
            mu <- exp(1.2+.6*sin(6*d$x)+b[d$g])
            d$y <- rnbinom(nrow(d),mu=mu,size=phi)
            args$family <- nbinom2()
            if (case == "zinb") {
                d$y[runif(nrow(d))<plogis(-.6+.8*cos(6*d$x))] <- 0
                args$ziformula <- ~s(x,k=5)
                args$dispformula <- ~s(x,k=5)
            } else args$dispformula <- ~s(z,k=5)
        }
        args$data <- d
        setup <- do.call(glmmTMB,c(args,list(doFit=FALSE)))
        old_start <- list(theta=rep(0,length(setup$parameters$theta)))
        a <- do.call(glmmTMB,args)
        old <- do.call(glmmTMB,c(args,list(start=old_start)))
        expect_equal(a$fit$convergence,0,info=case)
        expect_true(a$sdr$pdHess,info=case)
        expect_true(abs(a$fit$objective-old$fit$objective)<1e-4,info=case)
    }
})

test_that("inverse-link Gamma mean starts are not changed by smooth initialization", {
    args <- smooth_start_glm_data("Gamma_inverse")
    a <- do.call(glmmTMB,c(args,list(doFit=FALSE)))
    expect_true(all(a$parameters$beta == 1))
    expect_true(all(a$parameters$b == 1))
    ## These pre-existing mean starts can be invalid. Initializing smooth
    ## variances does not attempt to repair them or change other coefficients.
    eta <- drop(a$data.tmb$X %*% a$parameters$beta + a$data.tmb$Z %*% a$parameters$b)
    expect_true(any(eta<=0))
})
