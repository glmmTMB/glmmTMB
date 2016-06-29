## test more exotic familes/model types

stopifnot(require("testthat"),
          require("glmmTMB"))

context("fitting exotic families")
test_that("beta", {
     set.seed(101)
     beta <- c(2,1)
     sd.re <- 1
     phi <- 0.1
     ngrp <- 10
     nobs <- 200
     eps <- 0.001
     x <- rnorm(nobs)
     f <- factor(rep(1:ngrp,nobs/ngrp))
     u <- rnorm(ngrp,sd=sd.re)
     eta <- beta[1]+beta[2]*x+u[f]
     mu <- plogis(eta)
     y <- rbeta(nobs,shape1=mu/phi,shape2=(1-mu)/phi)
     y <- pmin(1-eps,pmax(eps,y))
     dd <<- data.frame(x,y,f) ## global assignment for testthat
     m1 <- glmmTMB(y~x+(1|f),family=list(family="beta",link="logit"),
                   data=dd)
     expect_equal(fixef(m1)[[1]],
                  structure(c(1.98250567574413, 0.843382531038295),
                            .Names = c("(Intercept)", "x")),
                  tol=1e-5)
     expect_equal(c(VarCorr(m1)[[1]][[1]]),
                  0.432365330831596, tol=1e-5)
 })

test_that("nbinom", {
     set.seed(101)
     beta <- c(2,1)
     sd.re <- 5
     phi <- 0.1
     ngrp <- 10
     nobs <- 200
     x <- rnorm(nobs)
     f <- factor(rep(1:ngrp,nobs/ngrp))
     u <- rnorm(ngrp,sd=sd.re)
     eta <- beta[1]+beta[2]*x+u[f]
     mu <- exp(eta)
     y <- rnbinom(nobs,size=phi,mu=mu)
     dd <<- data.frame(x,y,f) ## global assignment for testthat
     
     m1 <- glmmTMB(y~x+(1|f),family=list(family="nbinom2",link="log"),
                   data=dd)
     expect_equal(fixef(m1)[[1]],
                  structure(c(2.09866748794435, 1.12703589660625),
                            .Names = c("(Intercept)", "x")),
                  tol=1e-5)
     expect_equal(c(VarCorr(m1)[[1]][[1]]),
                  9.52772758529216, tol=1e-5)
     expect_equal(sigma(m1),0.09922738,tol=1e-5)

     ## nbinom1
     ## to simulate, back-calculate shape parameters for NB2 ...
     nbphi <- 2
     nbvar <- nbphi*mu
     ## V = mu*(1+mu/k) -> mu/k = V/mu-1 -> k = mu/(V/mu-1)
     k <- mu/(nbvar/mu - 1)
     y <- rnbinom(nobs,size=k,mu=mu)
     dd <<- data.frame(x,y,f) ## global assignment for testthat
     m1 <- glmmTMB(y~x+(1|f),family=list(family="nbinom1",link="log"),
                   data=dd)

     expect_equal(c(unname(c(fixef(m1)[[1]])),
                    c(VarCorr(m1)[[1]][[1]]),
                    sigma(m1)),
       c(1.93154240357181, 0.992776302432081,
         16.4150254955972, 1.00770603513152),
       tol=1e-5)

 })

test_that("truncated", {
    set.seed(101)
    z <- rnbinom(1000,size=2,mu=exp(2))
    z <- z[z>0]
    g1 <- glmmTMB(z~1,family=list(family="truncated_nbinom2",
                            link="log"),
            data=data.frame(z))
    expect_equal(c(unname(fixef(g1)[[1]]),sigma(g1)),
                 c(1.980207,1.892970),tol=1e-5)
})
