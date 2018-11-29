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

    ## Compare 2-column and prop/size specification
    dd <- data.frame(success=1:10, failure=11:20)
    dd$size <- rowSums(dd)
    dd$prop <- local( success / size, dd)
    mod4 <- glmmTMB(cbind(success,failure)~1,family=binomial,data=dd)
    mod5 <- glmmTMB(prop~1,weights=size,family=binomial,data=dd)
    expect_equal( logLik(mod4)     , logLik(mod5) )
    expect_equal( fixef(mod4)$cond , fixef(mod5)$cond )

    ## Now with extra weights
    dd$w <- 2
    mod6 <- glmmTMB(cbind(success,failure)~1,family=binomial,data=dd,weights=w)
    mod7 <- glmmTMB(prop~1,weights=size*w,family=binomial,data=dd)
    mod6.glm <- glm(cbind(success,failure)~1,family=binomial,data=dd,weights=w)
    mod7.glm <- glm(prop~1,weights=size*w,family=binomial,data=dd)
    expect_equal( logLik(mod6)[[1]]     , logLik(mod6.glm)[[1]] )
    expect_equal( logLik(mod7)[[1]]     , logLik(mod7.glm)[[1]] )
    expect_equal( fixef(mod6)$cond , fixef(mod7)$cond )
		
    ## Test TRUE/FALSE specification
    x <- c(TRUE, TRUE, FALSE)
    m1 <- glmmTMB(x~1, family=binomial())
    m2 <- glm    (x~1, family=binomial())
    expect_equal(
        as.numeric(logLik(m1)),
        as.numeric(logLik(m2))
    )
    expect_equal(
        as.numeric(unlist(fixef(m1))),
        as.numeric(coef(m2))
    )

    ## Mis-specifications
    prop <- c(.1, .2, .3)  ## weights=1 => prop * weights non integers
    expect_warning( glmmTMB(prop~1, family=binomial()) )   ## Warning as glm
    x <- c(1, 2, 3)        ## weights=1 => x > weights !
    expect_error  ( glmmTMB(x~1, family=binomial()) )      ## Error as glm
})
context("fitting exotic families")
test_that("beta", {
    set.seed(101)
    nobs <- 200; eps <- 0.001; phi <- 0.1
    dd0 <- simfun0(nobs=nobs,sd.re=1,invlink=plogis)
    y <- with(dd0,rbeta(nobs,shape1=mu/phi,shape2=(1-mu)/phi))
    dd <<- data.frame(dd0,y=pmin(1-eps,pmax(eps,y)))
    m1 <- glmmTMB(y~x+(1|f),family=beta_family(),
                  data=dd)
    expect_equal(fixef(m1)[[1]],
                 structure(c(1.98250567574413, 0.843382531038295),
                           .Names = c("(Intercept)", "x")),
                 tol=1e-5)
    expect_equal(c(VarCorr(m1)[[1]][[1]]),
                 0.433230926800709, tol=1e-5)
    ## allow family="beta", but with warning
    expect_warning(m2 <- glmmTMB(y~x+(1|f),family="beta",
                  data=dd),"please use")
    expect_equal(coef(summary(m1)),coef(summary(m2)))
    
 })

test_that("nbinom", {
    nobs <- 200; phi <- 0.1
    set.seed(101)
    dd0 <- simfun0(nobs=nobs)
    ## global assignment for testthat (??)
    dd <- data.frame(dd0,y=rnbinom(nobs,size=phi,mu=dd0$mu))
    m1 <- glmmTMB(y~x+(1|f),family=nbinom2(),
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
     m1 <- glmmTMB(y~x+(1|f),family=nbinom1(),
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
    expect_equal(fixef(glmmTMB(obs~x, family=nbinom2(link="identity"), dat)),
       structure(list(cond = structure(c(115.092240041138, 1.74390840106971),
       .Names = c("(Intercept)", "x")), zi = numeric(0),
       disp = structure(1.71242627201796, .Names = "(Intercept)")),
       .Names = c("cond", "zi", "disp"), class = "fixef.glmmTMB"))


    ## segfault (GH #248)
    dd <- data.frame(success=1:10,failure=10)
    expect_error(glmmTMB(cbind(success,failure)~1,family=nbinom2,data=dd),
                 "matrix-valued responses are not allowed")
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
                  family=betabinomial(),
                  data=dd)
    expect_equal(c(unname(c(fixef(m1)[[1]])),
                   c(VarCorr(m1)[[1]][[1]]),
                   sigma(m1)),
                 c(2.1482114,1.0574946,0.7016553,8.3768711),
                 tolerance=1e-5)
    ## Two-column specification
    m2 <- glmmTMB(cbind(y, N-y) ~ x + (1|f),
                  family=betabinomial(),
                  data=dd)
    expect_identical(m1$fit, m2$fit)
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
    g1_tp <- glmmTMB(z_tp~1,family=truncated_poisson(),
                  data=data.frame(z_tp))
    expect_equal(unname(fixef(g1_tp)[[1]]),0.9778593,tol=1e-5)
    ## Truncated poisson with zeros => invalid:
    num_zeros <- 10
    z_tp0 <<- c(rep(0, num_zeros), z_tp)
    expect_error(g1_tp0 <- glmmTMB(z_tp0~1,family=truncated_poisson(),
                      data=data.frame(z_tp0)))
    ## Truncated poisson with zeros and zero-inflation:
    g1_tp0 <- glmmTMB(z_tp0~1,family=truncated_poisson(),
                      ziformula=~1,
                      data=data.frame(z_tp0))
    expect_equal( plogis(as.numeric(fixef(g1_tp0)$zi)), num_zeros/length(z_tp0), tol=1e-7 ) ## Test zero-prob
    expect_equal(fixef(g1_tp0)$cond,  fixef(g1_tp)$cond, tol=1e-6) ## Test conditional model
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
    g1_nb2 <- glmmTMB(z_nb~1,family=truncated_nbinom2(),
            data=data.frame(z_nb))
    expect_equal(c(unname(fixef(g1_nb2)[[1]]),sigma(g1_nb2)),
                 c(1.980207,1.892970),tol=1e-5)
    ## Truncated nbinom2 with zeros => invalid:
    num_zeros <- 10
    z_nb0 <<- c(rep(0, num_zeros), z_nb)
    expect_error(g1_nb0 <- glmmTMB(z_nb0~1,family=truncated_nbinom2(),
                      data=data.frame(z_nb0)))
    ## Truncated nbinom2 with zeros and zero-inflation:
    g1_nb0 <- glmmTMB(z_nb0~1,family=truncated_nbinom2(),
                      ziformula=~1,
                      data=data.frame(z_nb0))
    expect_equal( plogis(as.numeric(fixef(g1_nb0)$zi)), num_zeros/length(z_nb0), tol=1e-7 ) ## Test zero-prob
    expect_equal(fixef(g1_nb0)$cond, fixef(g1_nb2)$cond, tol=1e-6) ## Test conditional model
    ## nbinom1: constant mean, so just a reparameterization of
    ##     nbinom2 (should have the same likelihood)
    ## phi=(1+mu/k)=1+exp(2)/2 = 4.69
    if (FALSE) {
        ## library(glmmADMB)
        g0_nb1 <- glmmadmb(z_nb~1,family="truncnbinom1",link="log")
        fixef(g0_nb1) ## 2.00112
        g0_nb1$alpha ## 3.784
    }
    g1_nb1 <- glmmTMB(z_nb~1,family=truncated_nbinom1(),
            data=data.frame(z_nb))
    expect_equal(c(unname(fixef(g1_nb1)[[1]]),sigma(g1_nb1)),
                 c(1.980207,3.826909),tol=1e-5)
    ## Truncated nbinom1 with zeros => invalid:
    expect_error(g1_nb0 <- glmmTMB(z_nb0~1,family=truncated_nbinom1(),
                      data=data.frame(z_nb0)))
    ## Truncated nbinom2 with zeros and zero-inflation:
    g1_nb0 <- glmmTMB(z_nb0~1,family=truncated_nbinom1(),
                      ziformula=~1,
                      data=data.frame(z_nb0))
    expect_equal( plogis(as.numeric(fixef(g1_nb0)$zi)), num_zeros/length(z_nb0), tol=1e-7 ) ## Test zero-prob
    expect_equal(fixef(g1_nb0)$cond, fixef(g1_nb1)$cond, tol=1e-6) ## Test conditional model
})

##Genpois
test_that("truncated_genpois",{
    tgp1 <<- glmmTMB(z_nb ~1, data=data.frame(z_nb), family=truncated_genpois())
    tgpdat <<- data.frame(y=simulate(tgp1)[,1])
    tgp2 <<- glmmTMB(y ~1, tgpdat, family=truncated_genpois())
    expect_equal(sigma(tgp1), sigma(tgp2), tol=1e-1)
    expect_equal(fixef(tgp1)$cond[1], fixef(tgp2)$cond[1], tol=1e-2)
    expect_lt(confint(tgp2)["sigma", "2.5 %"], sigma(tgp1))
    expect_lt(sigma(tgp1), confint(tgp2)["sigma", "97.5 %"])
    expect_lt(confint(tgp2)["cond.(Intercept)", "2.5 %"], unname(fixef(tgp1)$cond[1]))
    expect_lt(unname(fixef(tgp1)$cond[1]), confint(tgp2)["cond.(Intercept)", "97.5 %"])
})


context("trunc compois")
##Compois
test_that("truncated_compois",{
	cmpdat <<- data.frame(f=factor(rep(c('a','b'), 10)),
	 			y=c(15,5,20,7,19,7,19,7,19,6,19,10,20,8,21,8,22,7,20,8))
	tcmp1 <<- glmmTMB(y~f, cmpdat, family= truncated_compois())
	expect_equal(unname(fixef(tcmp1)$cond), c(2.9652730653, -0.9773987194), tol=1e-6)
	expect_equal(sigma(tcmp1), 0.1833339, tol=1e-6)
	expect_equal(predict(tcmp1,type="response")[1:2], c(19.4, 7.3), tol=1e-6)    
})

context("compois")
test_that("compois", {
#	cmpdat <<- data.frame(f=factor(rep(c('a','b'), 10)),
#	 			y=c(15,5,20,7,19,7,19,7,19,6,19,10,20,8,21,8,22,7,20,8))
	cmp1 <<- glmmTMB(y~f, cmpdat, family=compois())
	expect_equal(unname(fixef(cmp1)$cond), c(2.9652730653, -0.9773987194), tol=1e-6)
	expect_equal(sigma(cmp1), 0.1833339, tol=1e-6)
	expect_equal(predict(cmp1,type="response")[1:2], c(19.4, 7.3), tol=1e-6)
})

context("genpois")
test_that("genpois", {
	gendat <<- data.frame(y=c(11,10,9,10,9,8,11,7,9,9,9,8,11,10,11,9,10,7,13,9))
	gen1 <<- glmmTMB(y~1, family=genpois(), gendat)
	expect_equal(unname(fixef(gen1)$cond), 2.251292, tol=1e-6)
	expect_equal(sigma(gen1), 0.235309, tol=1e-6)
})

context("tweedie")
test_that("tweedie", {
    ## Boiled down tweedie:::rtweedie :
    rtweedie <- function (n, xi = power, mu, phi, power = NULL)
    {
        mu <- array(dim = n, mu)
        if ((power > 1) & (power < 2)) {
            rt <- array(dim = n, NA)
            lambda <- mu^(2 - power)/(phi * (2 - power))
            alpha <- (2 - power)/(1 - power)
            gam <- phi * (power - 1) * mu^(power - 1)
            N <- rpois(n, lambda = lambda)
            for (i in (1:n)) {
                rt[i] <- sum(rgamma(N[i], shape = -alpha, scale = gam[i]))
            }
        } else stop()
        as.vector(rt)
    }
    ## Simulation experiment
    nobs <- 2000; mu <- 4; phi <- 2; p <- 1.7
    set.seed(101)
    y <- rtweedie(nobs, mu=mu, phi=phi, power=p)
    twm <- glmmTMB(y ~ 1, family=tweedie())
    ## Check mu
    expect_equal(unname( exp(fixef(twm)$cond) ),
                 mu,
                 tolerance = .1)
    ## Check phi
    expect_equal(unname( exp(fixef(twm)$disp) ),
                 phi,
                 tolerance = .1)
    ## Check power
    expect_equal(unname( plogis(twm$fit$par["thetaf"]) + 1 ),
                 p,
                 tolerance = .01)
    ## Check internal rtweedie used by simulate
    y2 <- c(simulate(twm)[,1],simulate(twm)[,1])
    twm2 <- glmmTMB(y2 ~ 1, family=tweedie())
    expect_equal(fixef(twm)$cond, fixef(twm2)$cond, tol=1e-1)
    expect_equal(sigma(twm), sigma(twm2), tol=1e-1)
})

context("link function info available")

fam1 <- c("poisson","nbinom1","nbinom2","compois")
fam2 <- c("binomial","beta_family","betabinomial","tweedie")
for (f in c(fam1,paste0("truncated_",fam1),fam2)) {
    ## print(f)
    expect_true("linkinv" %in% names(get(f)()))
}

context("link info added to family")

d.AD <- data.frame(counts=c(18,17,15,20,10,20,25,13,12),
                   outcome=gl(3,1,9),
                   treatment=gl(3,3))
glm.D93 <- glmmTMB(counts ~ outcome + treatment, family = poisson(),
                   d.AD)
expect_warning(glm.D93B <- glmmTMB(counts ~ outcome + treatment,
                    family = list(family="poisson", link="log"),
                    d.AD))
## note update(..., family= ...) is only equal up to tolerance=5e-5 ...
glm.D93C <- glmmTMB(counts ~ outcome + treatment,
                    family = "poisson",
                    d.AD)
expect_equal(predict(glm.D93),predict(glm.D93B))
expect_equal(predict(glm.D93),predict(glm.D93C))

