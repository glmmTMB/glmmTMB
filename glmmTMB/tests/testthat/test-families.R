## test more exotic familes/model types

stopifnot(require("testthat"),
          require("glmmTMB"))

simfun0 <- function(beta=c(2,1),
                   sd.re=5,
                   ngrp=10,nobs=200,
                   invlink=exp) {
    x <- rnorm(nobs)
    f <- factor(rep(1:ngrp,nobs/ngrp))
    u <- rnorm(ngrp,sd=sd.re)
    eta <- beta[1]+beta[2]*x+u[f]
    mu <- invlink(eta)
    return(data.frame(x,f,mu))
}

context("alternative binomial specifications")
test_that("binomial", {
    load(system.file("testdata","radinger_dat.RData",package="lme4"))
    radinger_dat <<- radinger_dat ## global assignment for testthat
    mod1 <<- glmmTMB(presabs~predictor+(1|species),family=binomial,
                    radinger_dat)
    mod2 <<- update(mod1,as.logical(presabs)~.)
    expect_equal(predict(mod1),predict(mod2))
})
context("fitting exotic families")
test_that("beta", {
    set.seed(101)
    nobs <- 200; eps <- 0.001; phi <- 0.1
    dd0 <- simfun0(nobs=nobs,sd.re=1,invlink=plogis)
    y <- with(dd0,rbeta(nobs,shape1=mu/phi,shape2=(1-mu)/phi))
    dd <<- data.frame(dd0,y=pmin(1-eps,pmax(eps,y)))
    m1 <- glmmTMB(y~x+(1|f),family=list(family="beta",link="logit"),
                  data=dd)
    expect_equal(fixef(m1)[[1]],
                 structure(c(1.98250567574413, 0.843382531038295),
                           .Names = c("(Intercept)", "x")),
                 tol=1e-5)
    expect_equal(c(VarCorr(m1)[[1]][[1]]),
                 0.433230926800709, tol=1e-5)
 })

test_that("nbinom", {
    nobs <- 200; phi <- 0.1
    set.seed(101)
    dd0 <- simfun0(nobs=nobs)
    ## global assignment for testthat (??)
    dd <- data.frame(dd0,y=rnbinom(nobs,size=phi,mu=dd0$mu))
    m1 <- glmmTMB(y~x+(1|f),family=list(family="nbinom2",link="log"),
                  data=dd)
    expect_equal(fixef(m1)[[1]],
                 structure(c(2.09866748794435, 1.12703589660625),
                           .Names = c("(Intercept)", "x")),
                 tol=1e-5)
     expect_equal(c(VarCorr(m1)[[1]][[1]]),
                  9.54680210862774, tol=1e-5)
     expect_equal(sigma(m1),0.09922738,tol=1e-5)

     ## nbinom1
     ## to simulate, back-calculate shape parameters for NB2 ...
     nbphi <- 2
     nbvar <- nbphi*dd0$mu  ## n.b. actual model is (1+phi)*var,
                        ## so estimate of phi is approx. 1
     ## V = mu*(1+mu/k) -> mu/k = V/mu-1 -> k = mu/(V/mu-1)
     k <- with(dd0,mu/(nbvar/mu - 1))
     y <- rnbinom(nobs,size=k,mu=dd$mu)
     dd <- data.frame(dd0,y=y) ## global assignment for testthat
     m1 <- glmmTMB(y~x+(1|f),family=list(family="nbinom1",link="log"),
                   data=dd)
     expect_equal(c(unname(c(fixef(m1)[[1]])),
                    c(VarCorr(m1)[[1]][[1]]),
                    sigma(m1)),
       c(1.93154240357181, 0.992776302432081,
         16.447888398429, 1.00770603513152),
       tol=1e-5)

    ## identity link: GH #20
    x <- 1:100; m <- 2; b <- 100
    y <- m*x+b
    set.seed(101)
    dat <<- data.frame(obs=rnbinom(length(y), mu=y, size=5), x=x)
    ## with(dat, plot(x, obs))
    ## coef(mod1 <- MASS::glm.nb(obs~x,link="identity",dat))
    expect_equal(fixef(glmmTMB(obs~x, family=list(family="nbinom2",
                                             link="identity"), dat)),
       structure(list(cond = structure(c(115.092240041138, 1.74390840106971),
       .Names = c("(Intercept)", "x")), zi = numeric(0),
       disp = structure(1.71242627201796, .Names = "(Intercept)")),
       .Names = c("cond", "zi", "disp"), class = "fixef.glmmTMB"))

 })

test_that("dbetabinom", {
    set.seed(101)
    nobs <- 200; eps <- 0.001; phi <- 0.1
    dd0 <- simfun0(nobs=nobs,sd.re=1,invlink=plogis)
    p <- with(dd0,rbeta(nobs,shape1=mu/phi,shape2=(1-mu)/phi))
    p <- pmin(1-eps,pmax(p,eps))
    b <- rbinom(nobs,size=5,prob=p)
    dd <<- data.frame(dd0,y=b,N=5)
    m1 <- glmmTMB(y/N~x+(1|f),
                  weights=N,
                  family=list(family="betabinomial",link="logit"),
                  data=dd)
    expect_equal(c(unname(c(fixef(m1)[[1]])),
                   c(VarCorr(m1)[[1]][[1]]),
                   sigma(m1)),
                 c(2.1482114,1.0574946,0.7016553,8.3768711),
                 tolerance=1e-5)
})

test_that("truncated", {
    ## Poisson
    set.seed(101)
    z_tp <<- rpois(1000,lambda=exp(1))
    z_tp <<- z_tp[z_tp>0]
    if (FALSE) {
        ## n.b.: keep library() calls commented out, they may
        ##   trigger CRAN complaints
        ## library(glmmADMB)
        g0_tp <- glmmadmb(z_tp~1,family="truncpoiss",link="log")
        fixef(g0) ## 0.9778591
    }
    g1_tp <- glmmTMB(z_tp~1,family=list(family="truncated_poisson",
                            link="log"),
                  data=data.frame(z_tp))
    expect_equal(unname(fixef(g1_tp)[[1]]),0.9778593,tol=1e-5)
    ## nbinom2
    set.seed(101)
    z_nb <<- rnbinom(1000,size=2,mu=exp(2))
    z_nb <<- z_nb[z_nb>0]
    if (FALSE) {
        ## library(glmmADMB)
        g0_nb2 <- glmmadmb(z_nb~1,family="truncnbinom",link="log")
        fixef(g0_nb2) ## 1.980207
        g0_nb2$alpha ## 1.893
    }
    g1_nb2 <- glmmTMB(z_nb~1,family=list(family="truncated_nbinom2",
                            link="log"),
            data=data.frame(z_nb))
    expect_equal(c(unname(fixef(g1_nb2)[[1]]),sigma(g1_nb2)),
                 c(1.980207,1.892970),tol=1e-5)
    ## nbinom1: constant mean, so just a reparameterization of
    ##     nbinom2 (should have the same likelihood)
    ## phi=(1+mu/k)=1+exp(2)/2 = 4.69
    if (FALSE) {
        ## library(glmmADMB)
        g0_nb1 <- glmmadmb(z_nb~1,family="truncnbinom1",link="log")
        fixef(g0_nb1) ## 2.00112
        g0_nb1$alpha ## 3.784
    }
    g1_nb1 <- glmmTMB(z_nb~1,family=list(family="truncated_nbinom1",
                            link="log"),
            data=data.frame(z_nb))
    expect_equal(c(unname(fixef(g1_nb1)[[1]]),sigma(g1_nb1)),
                 c(1.980207,3.826909),tol=1e-5)
})
test_that("compois", {
	cmpdat <<- data.frame(f=factor(rep(c('a','b'), 10)),
	 			y=c(15,5,20,7,19,7,19,7,19,6,19,10,20,8,21,8,22,7,20,8))
	cmp1 <<- glmmTMB(y~f, cmpdat, family="compois")
	expect_equal(unname(fixef(cmp1)$cond), c(2.9652731, -0.9773987), tol=1e-6)
	expect_equal(sigma(cmp1), 0.1833339, tol=1e-6)
	expect_equal(predict(cmp1)[1:2], c(19.4, 7.3), tol=1e-6)
})
test_that("genpois", {
	gendat <<- data.frame(y=c(11,10,9,10,9,8,11,7,9,9,9,8,11,10,11,9,10,7,13,9))
	gen1 <<- glmmTMB(y~1, family="genpois", gendat)
	expect_equal(unname(fixef(gen1)$cond), 2.251292, tol=1e-6)
	expect_equal(sigma(gen1), 0.235309, tol=1e-6)
})