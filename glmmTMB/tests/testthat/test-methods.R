stopifnot(require("testthat"),
          require("glmmTMB"))

data(sleepstudy, cbpp, Pastes,
     package = "lme4")

if (getRversion() < "3.3.0") {
    sigma.default <- function (object, use.fallback = TRUE, ...) 
        sqrt(deviance(object, ...)/(nobs(object, use.fallback = use.fallback) - 
                                    length(coef(object))))
}

load(system.file("test_data", "models.rda", package="glmmTMB",
                 mustWork=TRUE))

context("basic methods")

test_that("Fitted and residuals", {
    expect_equal(length(fitted(fm2)),nrow(sleepstudy))
    expect_equal(mean(fitted(fm2)),298.507891)
    expect_equal(mean(residuals(fm2)),0,tol=1e-5)
    ## Pearson and response are the same for a Gaussian model
    expect_equal(residuals(fm2,type="response"),
                 residuals(fm2,type="pearson"))
    ## ... but not for Poisson or NB ...
    expect_false(mean(residuals(fm2P,type="response"))==
                 mean(residuals(fm2P,type="pearson")))
    expect_false(mean(residuals(fm2NB,type="response"))==
                 mean(residuals(fm2NB,type="pearson")))
    rr2 <- function(x) sum(residuals(x,type="pearson")^2)
    ## test Pearson resids for gaussian, Gamma vs. base-R versions
    ss <- as.data.frame(state.x77)
    expect_equal(rr2(glm(Murder~Population,ss,family=gaussian)),
          rr2(glmmTMB(Murder~Population,ss,family=gaussian)))
    expect_equal(rr2(glm(Murder~Population,ss,family=Gamma(link="log"))),
                 rr2(glmmTMB(Murder~scale(Population),ss,
                             family=Gamma(link="log"))),tol=1e-5)
    ## weights are incorporated in Pearson residuals
    ## GH 307
    tmbm4 <- glm(incidence/size ~ period,
             data = cbpp, family = binomial, weights = size)
    tmbm5 <- glmmTMB(incidence/size ~ period,
                     data = cbpp, family = binomial, weights = size)
    expect_equal(residuals(tmbm4,type="pearson"),
                 residuals(tmbm5,type="pearson"),tolerance=1e-6)
    ## two-column responses give vector of residuals GH 307
    tmbm6 <- glmmTMB(cbind(incidence,size-incidence) ~ period,
                     data = cbpp, family = binomial)
    expect_equal(residuals(tmbm4,type="pearson"),
                 residuals(tmbm6,type="pearson"),tolerance=1e-6)

    ## predict handles na.exclude correctly
    ## GH 568
    b <- rnorm(332)
    mu <- exp(1.5 + .26*b)
    y <- sapply(mu, function(mu){rpois(1, lambda = mu)})
    napos <- 51
    b[napos] <- NA
    y.na <- y
    y.na[napos] <- NA 
    mod.ex <- glmmTMB(y ~ b, family = "poisson", na.action = "na.exclude")
    ## Get predictions/resids
    pr.ex <- predict(mod.ex, type = "response") # SEEMS to work fine
    expect_equal(which(is.na(pr.ex)),napos)
    rs.ex <- residuals(mod.ex, type = "response")
    expect_equal(unname(which(is.na(rs.ex))),napos)
    pr.rs.ex <- pr.ex + rs.ex
    expect_equal(unname(pr.rs.ex), y.na)
})

test_that("Predict", {
    expect_equal(predict(fm2),predict(fm2,newdata=sleepstudy))
    pr2se <- predict(fm2, se.fit=TRUE)
    i <- sample(nrow(sleepstudy), 20)
    newdata <- sleepstudy[i, ]
    pr2sub <- predict(fm2, newdata=newdata, se.fit=TRUE)
    expect_equivalent(pr2se$fit, predict(fm2))
    expect_equivalent(pr2se$fit[i], pr2sub$fit)
    expect_equivalent(pr2se$se.fit[i], pr2sub$se.fit)
    expect_equal(unname( pr2se$   fit[1] ), 254.2208, tol=1e-4)
    expect_equal(unname( pr2se$se.fit[1] ), 12.94514, tol=1e-4)
    expect_equal(unname( pr2se$   fit[100] ), 457.9684, tol=1e-4)
    expect_equal(unname( pr2se$se.fit[100] ), 14.13943, tol=1e-4)

    ## predict without response in newdata
    expect_equal(predict(fm2),
                 predict(fm2,newdata=sleepstudy[,c("Days","Subject")]))
    
})


test_that("VarCorr", {
   vv <- VarCorr(fm2)
   vv2 <- vv$cond$Subject
   expect_equal(dim(vv2),c(2,2))
   expect_equal(outer(attr(vv2,"stddev"),
                      attr(vv2,"stddev"))*attr(vv2,"correlation"),
                vv2,check.attributes=FALSE)
   vvd <- VarCorr(fm2diag)
   expect_equal(vvd$cond$Subject[1,2],0) ## off-diagonal==0
})

test_that("drop1", {
      dd <- drop1(fm2,test="Chisq")
      expect_equal(dd$AIC,c(1763.94,1785.48),tol=1e-4)              
          })
test_that("anova", {
      aa <- anova(fm0,fm2)
      expect_equal(aa$AIC,c(1785.48,1763.94),tol=1e-4)
          })

test_that("terms", {
    ## test whether terms() are returned with predvars for doing
    ## model prediction etc. with complex bases
    dd <<- data.frame(x=1:10,y=1:10)
    require("splines")
    ## suppress convergence warnings(we know this is a trivial example)
    suppressWarnings(m <- glmmTMB(y~ns(x,3),dd))
    ## if predvars is not properly attached to term, this will
    ## fail as it tries to construct a 3-knot spline from a single point
    expect_equal(model.matrix(delete.response(terms(m)),data=data.frame(x=1)),
      structure(c(1, 0, 0, 0), .Dim = c(1L, 4L), .Dimnames = list("1", 
    c("(Intercept)", "ns(x, 3)1", "ns(x, 3)2", "ns(x, 3)3")),
    assign = c(0L, 1L, 1L, 1L)))
})

test_that("terms back-compatibility", {
    skip("skip until TMB fits are back-compatible")
    f0 <- readRDS(system.file("test_data", "oldfit.rds",
                              package="glmmTMB",
                              mustWork=TRUE))
    expect_true(!is.null(terms(f0)))
})
    
test_that("summary_print", {
    getVal <- function(x,tag="Dispersion") {
        cc <- capture.output(print(summary(x)))
        if (length(gg <- grep(tag,cc,value=TRUE))==0) return(NULL)
        cval <- sub("^.*: ","",gg) ## get value after colon ...
        return(as.numeric(cval))
    }
    ## no dispersion printed for Gaussian or disp==1 families
    expect_equal(getVal(fm2),654.9,tolerance=1e-2)
    expect_equal(getVal(fm2P),NULL)
    expect_equal(getVal(fm2G),0.00654,tolerance=1e-2)
    expect_equal(getVal(fm2NB,"Dispersion"),286,tolerance=1e-2)
})

test_that("sigma", {
    s1 <<- sigma(lm(Reaction~Days,sleepstudy))
    s2 <<- sigma(glm(Reaction~Days,sleepstudy,family=Gamma(link="log")))
    s3 <<- MASS::glm.nb(round(Reaction)~Days,sleepstudy)
    ## remove bias-correction
    expect_equal(sigma(fm3),s1*(1-1/nobs(fm3)),tolerance=1e-3)
    expect_equal(sigma(fm3G),s2,tolerance=5e-3)
    expect_equal(s3$theta,sigma(fm3NB),tolerance=1e-4)
})

test_that("confint", {
    ci <- confint(fm2, 1:2, estimate=FALSE)
    expect_equal(ci,
        structure(c(238.406083254105, 7.52295734348693,
                    264.404107485727, 13.4116167530013),
                  .Dim = c(2L, 2L),
                  .Dimnames = list(c("(Intercept)", "Days"),
                                   c("2.5 %", "97.5 %"))),
        tolerance=1e-6)
    ciw <- confint(fm2, 1:2, method="Wald", estimate=FALSE)
    expect_warning(confint(fm2,type="junk"),
                   "extra arguments ignored")
    ## Gamma test Std.Dev and sigma
    ci.2G <- confint(fm2G, full=TRUE, estimate=FALSE)
    ci.2G.expect <- structure(c(5.4810173444768, 0.0247781468857994, 0.0676097043327788, 
                             0.0115949839191128, -0.518916570291726, 0.0720456818399729, 5.58401849115119, 
                             0.0429217639222305, 0.150456372618643, 0.0264376535768207, 0.481694558481224, 
                             0.0907365112123184),
                           .Dim = c(6L, 2L),
                           .Dimnames = list(c("cond.(Intercept)", 
                                              "cond.Days", "cond.Std.Dev.(Intercept)",
                                              "cond.Std.Dev.Days", "cond.Cor.Days.(Intercept)", "sigma"),
                                            c("2.5 %", "97.5 %")))
    expect_equal(ci.2G, ci.2G.expect, tolerance=1e-6)
    ## nbinom2 test Std.Dev and sigma
    ci.2NB <- confint(fm2NB, full=TRUE, estimate=FALSE)
    ci.2NB.expect <-
        structure(c(5.48098713986992, 0.0248163859092965, 0.066177247560203, 
                    0.0113436356932709, -0.520883841816814, 183.810584738707, 5.58422550782448, 
                    0.0428993227431795, 0.150917850214506, 0.026354988318893, 0.502211676507888, 
                    444.735668635694),
                  .Dim = c(6L, 2L),
                  .Dimnames = list(c("cond.(Intercept)", 
                                     "cond.Days",
                                     "cond.Std.Dev.(Intercept)", "cond.Std.Dev.Days", 
                                     "cond.Cor.Days.(Intercept)", "sigma"), c("2.5 %", "97.5 %")))
    expect_equal(ci.2NB, ci.2NB.expect, tolerance=1e-6)
    ## profile CI
    ## ... no RE
    ci.prof0 <- confint(fm_noRE, full=TRUE, method="profile", npts=3)
    expect_equal(ci.prof0,
                 structure(c(238.216039176535, 7.99674863649355, 7.51779308310198, 
                             264.368471102549, 12.8955469713508, 7.93347860201449),
                           .Dim = 3:2, .Dimnames = list(c("(Intercept)", "Days", "d~(Intercept)"),
                                                        c("2.5 %", "97.5 %"))),
                 tolerance=1e-5)

    ci.prof <- confint(fm2,parm=1,method="profile", npts=3)
    expect_equal(ci.prof,
                 structure(c(237.27249, 265.13383),
                           .Dim = 1:2, .Dimnames = list(
                                "(Intercept)", c("2.5 %", "97.5 %"))),
                 tolerance=1e-6)
    ## uniroot CI
    ci.uni <- confint(fm2,parm=1,method="uniroot")
    expect_equal(ci.uni,
                 structure(c(237.68071,265.12949,251.4050979),
                        .Dim = c(1L, 3L),
        .Dimnames = list("(Intercept)", c("2.5 %", "97.5 %", "Estimate"))),
                 tolerance=1e-6)
    ## check against 'raw' tmbroot
    tmbr <- TMB::tmbroot(fm2$obj,name=1)
    expect_equal(ci.uni[1:2],unname(c(tmbr)))

    ## GH #438
    cc <- confint(fm4)
    expect_equal(dim(cc),c(5,3))
    expect_equal(rownames(cc),
                 c("(Intercept)", "Illiteracy", "Population", "Area", "`HS Grad`"))

})

test_that("profile", {
    p1_th <- profile(fm1,parm="theta_",npts=4)
    expect_true(all(p1_th$.par=="theta_1|Subject.1"))
    p1_b <- profile(fm1,parm="beta_",npts=4)
    expect_equal(unique(as.character(p1_b$.par)),
                 c("(Intercept)","Days"))
})

test_that("profile (no RE)", {
    p0_th <- profile(fm_noRE,npts=4)
    expect_equal(dim(p0_th),c(43,3))
})

test_that("vcov", {
    expect_equal(dim(vcov(fm2)[[1]]),c(2,2))
    expect_equal(dim(vcov(fm2,full=TRUE)),c(6,6))
    expect_equal(rownames(vcov(fm2,full=TRUE)),
           structure(c("(Intercept)", "Days", "d~(Intercept)",
                       "theta_Days|Subject.1", "theta_Days|Subject.2",
                       "theta_Days|Subject.3"),
          .Names = c("cond1", "cond2", "disp", "theta1", "theta2", "theta3")))
    ## vcov doesn't include dispersion for non-dispersion families ...
    expect_equal(dim(vcov(fm2P,full=TRUE)),c(5,5))
    ## oops, dot_check() disabled in vcov.glmmTMB ...
    ## expect_error(vcov(fm2,x="junk"),"unknown arguments")
})

set.seed(101)
test_that("simulate", {
    sm2 <<- rowMeans(do.call(cbind, simulate(fm2, 10)))
    sm2P <<- rowMeans(do.call(cbind, simulate(fm2P, 10)))
    sm2G <<- rowMeans(do.call(cbind, simulate(fm2G, 10)))
    sm2NB <<- rowMeans(do.call(cbind, simulate(fm2NB, 10)))
    expect_equal(sm2, sleepstudy$Reaction, tol=20)
	expect_equal(sm2P, sleepstudy$Reaction, tol=20)
	expect_equal(sm2G, sleepstudy$Reaction, tol=20)
	expect_equal(sm2NB, sleepstudy$Reaction, tol=20)
})

test_that("formula", {
    expect_equal(formula(fm2),Reaction ~ Days + (Days | Subject))
    expect_equal(formula(fm2, fixed.only=TRUE),Reaction ~ Days)
    expect_equal(formula(fm2, component="disp"), ~1)
    expect_equal(formula(fm2, component="disp", fixed.only=TRUE), ~1)
    expect_equal(formula(fm2, component="zi"), ~0)
    expect_equal(formula(fm2, component="zi", fixed.only=TRUE), ~0)
})

context("simulate consistency with glm/lm")
test_that("binomial", {
    s1 <- simulate(f1b, 5, seed=1)
    s2 <- simulate(f2b, 5, seed=1)
    s3 <- simulate(f3b, 5, seed=1)
    expect_equal(max(abs(as.matrix(s1) - as.matrix(s2))), 0)
    expect_equal(max(abs(as.matrix(s1) - as.matrix(s3))), 0)
})

test_that("residuals from binomial factor responses", {
    expect_equal(residuals(fm2Bf),residuals(fm2Bn))
})

mkstr <- function(dd) {
    ff <- which(vapply(dd,is.factor,logical(1)))
    for (i in ff) {
        dd[[i]] <- as.character(dd[[i]])
    }
    return(dd)
}
rr <- function(txt) {
    read.table(header=TRUE,stringsAsFactors=FALSE,text=txt,
               colClasses=rep(c("character","numeric"),c(5,2)))
}
context("Ranef etc.")
test_that("as.data.frame(ranef(.)) works",
  {
      expect_equal(
          mkstr(as.data.frame(ranef(fm3ZIP))[c("cond.1","cond.19","zi.1"),]),
          rr(
"        component  grpvar        term grp       condval      condsd
cond.1       cond Subject (Intercept) 308  1.066599e-02 0.040430751
cond.19      cond Subject        Days 308  2.752424e-02 0.007036958
zi.1           zi Subject (Intercept) 308 -2.850238e-07 0.127106817
"),
tolerance=1e-5)
      expect_equal(
          mkstr(as.data.frame(ranef(fm2diag2))[c("cond.1","cond.19"),]),
          rr(
"        component  grpvar        term grp  condval    condsd
cond.1       cond Subject (Intercept) 308 1.854597 13.294388
cond.19      cond Subject        Days 308 9.236420  2.699692
"),
tolerance=1e-5)
  })

test_that("ranef(.) works with more than one grouping factor",
{
    expect_equal(sort(names(ranef(fmP)[["cond"]])), c("batch","sample"))
    expect_equal(dim(as.data.frame(ranef(fmP))), c(40,6))
})

test_that("coef(.) works", {
          cc <- coef(fm3ZIP)
          expect_equal(cc[["cond"]][[1]][1,],
                       structure(list(`(Intercept)` = 5.54291514202372,
                                      Days = 0.0613847280572168),
                                 row.names = "308", class = "data.frame"),
                       tolerance=1e-5)
          expect_equal(cc[["zi"]][[1]][1,,drop=FALSE],
                       structure(list(`(Intercept)` = -13.2478200379555), row.names = "308", class = "data.frame"),
                       tolerance=1e-5)
})

test_that("simplified coef(.) printing", {
    op <- options(digits=2)
    cc <- capture.output(print(coef(fm0)))          
    expect_equal(cc[1:3],c("$Subject", "    Days (Intercept)", "308 20.6         249"))
    options(op)
})

context("refit")

## weird stuff here with environments, testing ...
test_that("various binomial response types work", {
    ## FIXME: test for factors, explicit cbind(.,.)
    ## do we need to define this within this scope?
    ddb <- data.frame(y=I(yb))
    ddb <- within(ddb, {
        w <- rowSums(yb)
        prop <- y[,1]/w
    })
    s1 <- simulate(f1b, 1, seed=1)
    f1 <- fixef(refit(f1b,s1[[1]]))
    s3 <- simulate(f3b, 1, seed=1)
    f3 <- fixef(refit(f3b,s3[[1]]))
    expect_equal(f1,f3)
    expect_error(refit(f4b,s3[[1]]),
                  "can't find response in data")
})

test_that("binomial response types work with data in external scope", {
    s1 <- simulate(f1b, 1, seed=1)
    f1 <- fixef(refit(f1b,s1[[1]]))
    s3 <- simulate(f3b, 1, seed=1)
    f3 <- fixef(refit(f3b,s3[[1]]))
    expect_equal(f1,f3)
})

test_that("confint works for models with dispformula", {
    ## FIXME: should make this an example
    sim1 <- function(nfac=40, nt=100, facsd=0.1, tsd=0.15, mu=0, residsd=1) {
        dat <- expand.grid(fac=factor(letters[1:nfac]), t=1:nt)
        n <- nrow(dat)
        dat$REfac <- rnorm(nfac, sd=facsd)[dat$fac]
        dat$REt <- rnorm(nt, sd=tsd)[dat$t]
        dat$x <- rnorm(n, mean=mu, sd=residsd) + dat$REfac + dat$REt
        dat
    }
    set.seed(101)
    d1 <- sim1(mu=100, residsd=10)
    d2 <- sim1(mu=200, residsd=5)
    d1$sd <- "ten"
    d2$sd <- "five"
    dat <- rbind(d1, d2)
    m1 <- glmmTMB(x ~ sd + (1|t), dispformula=~sd, data=dat)
    ref_val <- structure(c(3.14851028784965, 1.30959944530366, 3.25722952319077, 
                           1.46335165911997, 3.20286990552021, 1.38647555221182), .Dim = 2:3,
                         .Dimnames = list(c("disp.(Intercept)", "disp.sdten"),
                                          c("2.5 %", "97.5 %", "Estimate")))
    expect_equal(tail(confint(m1),2), ref_val)
})

## utility functions for checking truncated-distribution simulations
## FIXME: add to utils.R?
dtruncated_nbinom2 <- function(x,size,mu,k=0,log=FALSE) {
    y <- ifelse(x<=k,-Inf,
                dnbinom(x,mu=mu, size=size,log=TRUE) -
                pnbinom(k, mu=mu, size=size, lower.tail=FALSE,
                        log.p=TRUE))
    if (log) return(y) else return(exp(y))
}

dtruncated_poisson <- function(x,lambda,k=0,log=FALSE) {
    y <- ifelse(x<=k,-Inf,
                dpois(x,lambda,log=TRUE) -
                ppois(k, lambda=lambda, lower.tail=FALSE,
                      log.p=TRUE))
    if (log) return(y) else return(exp(y))
}

dtruncated_nbinom1 <- function(x,phi,mu,k=0,log=FALSE) {
    ## V=mu*(1+phi) = mu*(1+mu/k) -> k=mu/phi
    size <- mu/phi
    y <- ifelse(x<=k,-Inf,
                dnbinom(x,mu=mu, size=size,log=TRUE) -
                pnbinom(k, mu=mu, size=size, lower.tail=FALSE,
                        log.p=TRUE))
    if (log) return(y) else return(exp(y))
}

simfun <- function(formula, family, data, beta=c(0,1)) {
    ss <- list(beta=beta)
    if (grepl("nbinom",family)) ss$betad <- 0
    suppressWarnings(m1 <- glmmTMB(formula,
                                   family=family,
                                   data=data,
                                   start=ss,
                                   control=glmmTMBControl(optCtrl=list(eval.max=0,iter.max=0))))
    return(m1)
}

ntab <- function(formula=y~x, family, data, seed=101) {
    set.seed(seed)
    m1 <- simfun(formula, family, data)
    return(table(exp(data$x),unlist(simulate(m1))))
}

pfun <- function(i,tab, dist="nbinom2", data) {
    n <- as.numeric(names(tab[i,]))
    plot(n,tab[i,]/sum(tab[i,]))
    m <- exp(data$x)[i]
    argList <- switch(dist,
                      nbinom1=list(n, phi=1, mu=m),
                      nbinom2=list(n, size=1, mu=m),
                      poisson=list(n, lambda=m))
    expected <- do.call(paste0("dtruncated_",dist), argList)
    lines(n,expected)
}

test_that("trunc nbinom simulation", {
    ## GH 572
    dd <- data.frame(f=factor(1:2),
                     y=rep(1,2))
    ## results for second element of sim, depending on family:
    simres <- c(truncated_nbinom2=1,truncated_nbinom1=2)
    for (f in paste0("truncated_nbinom",1:2)) {
        ## generate a model with two groups, one with a ridiculously low (log mean).
        ## don't allow the optimizer to actually do anything, so coefs will remain
        ## at their starting values
        m1 <- simfun(y~f, family=f, data=dd, beta=c(-40,39))
        expect_equal(fixef(m1)$cond, c(`(Intercept)` = -40, f2 = 39))
        expect_equal(fitted(m1),c(4.24835425529159e-18, 0.367879441171442))
        ## should NOT get NaN (or zero) for the first group if hack/fix is working
        expect_equal(unname(unlist(simulate(m1,seed=101))),c(1,1))
    }
    
})

test_that("trunc nbinom sim 2", {
    dd <- expand.grid(x=log(1:5),
                      rep=1:10000,
                      y=1)
    t1 <- ntab(family="truncated_nbinom1", data=dd)
    t2 <- ntab(family="truncated_nbinom2", data=dd)
    if (FALSE) { 
        op <- par(ask=TRUE)
        for (i in 1:nrow(t1)) pfun(i,tab=t1,dist="nbinom1",data=dd)
        for (i in 1:nrow(t2)) pfun(i,tab=t2,dist="nbinom2",data=dd)
        par(op)
    }

})

test_that("trunc poisson simulation", {
    dd <- expand.grid(x=log(1:5),
                      rep=1:10000,
                      y=1)
    t3 <- ntab(family="truncated_poisson", data=dd)
    expect_equal(unname(t3[1,1:6]),
                 c(5829L, 2905L, 963L, 242L, 56L, 5L))
    ## explore
    if (FALSE) { 
        op <- par(ask=TRUE)
        for (i in 1:nrow(t3)) pfun(i,tab=t3,dist="poisson",data=dd)
        par(op)
    }
})
