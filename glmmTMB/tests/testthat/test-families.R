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
