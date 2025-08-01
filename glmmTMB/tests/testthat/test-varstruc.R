stopifnot(require("testthat"),
          require("glmmTMB"),
          require("lme4"))
source(system.file("test_data/glmmTMB-test-funs.R",
                   package="glmmTMB", mustWork=TRUE))

data(sleepstudy, cbpp, package = "lme4")

## fsleepstudy and fitted models come from inst/test_data/models.rda
## OR running inst/test_data/make_ex.R

test_that("diag", {
   ## two formulations of diag and lme4 all give same log-lik
   expect_equal(logLik(fm_diag1),logLik(fm_diag2_lmer))
   expect_equal(logLik(fm_diag1),logLik(fm_diag2))
})

test_that("cs_us", {
## for a two-level factor, compound symmetry and unstructured
##  give same result
    expect_equal(logLik(fm_us1),logLik(fm_cs1))
    expect_equal(logLik(fm_us1),logLik(fm_us1_lmer))
})

test_that("cs_homog", {
    fm_csh   <- glmmTMB(Reaction ~ Days + homcs(Days | Subject), sleepstudy)
    vv <- unname(diag(VarCorr(fm_csh)[["cond"]][["Subject"]]))
    expect_identical(vv[1], vv[2])
})

test_that("cs_homog", {
## *homogenous* compound symmetry vs. nested random effects
    expect_equal(logLik(fm_nest),logLik(fm_nest_lmer))

})

test_that("basic ar1", {
    ## base fm_ar1 does include corr matrix
    vv <- VarCorr(fm_ar1)[["cond"]]
    expect_equal(attr(vv[[2]], "correlation")[2,1], 0.87299, tolerance = 1e-5)
    ## also need to test 
})

## change to something better behaved
test_that("print ar1 (>1 RE)", {
    ## sim order of sampling rnorm() values changed with implementation of hetar1, so use stored sim
    ##
    ## fsleepstudy$sim <- simulate_new(~ 1 + (1|Subject) + ar1(row+0| Subject),
    ##                                 newdata=fsleepstudy,
    ##                                 newparams = list(beta=0, betadisp = 1, theta = c(1, 1, 1)),
    ##                                 family = gaussian,
    ##                                 seed = 101)[[1]]
    ##     saveRDS(fsleepstudy$sim, file="../../inst/test_data/sim_ar1.rds",version=2)
    fsleepstudy$sim <- readRDS(system.file("test_data", "sim_ar1.rds", package="glmmTMB"))
    fm_ar2 <- glmmTMB(sim ~ 1 +
                          (1|Subject) + ar1(row+0| Subject), fsleepstudy)
    cco <- gsub(" +"," ",
                trimws(capture.output(print(summary(fm_ar2),digits=2))))
    expect_equal(cco[12:14],
                 c("Subject (Intercept) 7.0 2.6", "Subject.1 row1 5.9 2.4 0.78 (ar1)", 
                   "Residual 8.1 2.8"))
})

test_that("ar1 and hetar1 require factor time", {
  skip_on_cran()
    expect_error(glmmTMB(Reaction ~ 1 +
                             (1|Subject) + ar1(as.numeric(row)+0| Subject), fsleepstudy),
                 "expects a single")
    expect_error(glmmTMB(Reaction ~ 1 +
                             (1|Subject) + hetar1(as.numeric(row)+0| Subject), fsleepstudy),
                 "expects a single")
    ## works even when the factor is a weird/hard-to-recognize component
    expect_is(glmmTMB(Reaction ~ 1 +
                          (1|Subject) + ar1(relevel(factor(row),"2")+0| Subject),
                      fsleepstudy),
              "glmmTMB")
    expect_is(glmmTMB(Reaction ~ 1 +
                          (1|Subject) + hetar1(relevel(fDays,"[6,10)")+0| Subject),
                      fsleepstudy),
              "glmmTMB")
})

## FIXME: simpler to check formatVC() directly?
get_vcout <- function(x, g="\\bSubject\\b") {
    cc <- capture.output(print(VarCorr(x)))
    cc1 <- grep(g, cc, value=TRUE, perl=TRUE)
    ss <- strsplit(cc1,"[^[:alnum:][:punct:]]+")[[1]]
    return(ss[nchar(ss)>0])
}

test_that("varcorr_print", {
    skip_on_cran()
    ss <- get_vcout(fm_cs1)
    expect_equal(length(ss),5)
    expect_equal(ss[4:5],c("0.081","(cs)"))
    ss2 <- get_vcout(fm_ar1,g="\\bSubject.1\\b")
    expect_equal(length(ss2),5)
    expect_equal(ss2[4:5],c("0.873","(ar1)"))

    ## test case with two different size V-C
    set.seed(101)
    dd <- data.frame(y=rnorm(1000),c=factor(rep(1:2,500)),
                 w=factor(rep(1:10,each=100)),
                 s=factor(rep(1:10,100)))
    ## non-pos-def case (we don't care at the moment)
    m1 <- suppressWarnings(glmmTMB(y~c+(c|w)+(1|s),data=dd,
                  family=gaussian))
    cc <- squash_white(capture.output(print(VarCorr(m1),digits=2)))
    ## updated for var -> SD reparam
    expect_equal(cc,
                 c("Conditional model:", "Groups Name Std.Dev. Corr",
                   "w (Intercept) 9.6e-05", 
                   "c2 4.0e-06 0.99", "s (Intercept) 9.4e-06",
                   "Residual 9.6e-01"))

    ## check that all std devs are being printed (GH #851)
    cc <- capture.output(VarCorr(fm_cs2))
    expect_equal(length(cc), 7)
    expect_equal(length(grep("fDays", cc)), 2)
})

ff <- system.file("test_data","cov_struct_order.rds",package="glmmTMB")

if (nchar(ff)>0) {
    dat <- readRDS(ff)
} else {
    set.seed(101)
    nb <- 100
    ns <- nb*3
    nt <- 100
    cor <- .7
    dat  <-  data.frame(Block = factor(rep(1:nb, each = ns/nb*nt)),
                        Stand = factor(rep(1:ns, each = nt)),
                        Time = rep(1:nt, times = ns),
                        blockeff = rep(rnorm(nb, 0, .5), each = ns/nb*nt),
                        standeff = rep(rnorm(ns, 0, .8), each = nt),
                        resid = c(t(MASS::mvrnorm(ns, mu = rep(0, nt),
                                                  Sigma = 1.2*cor^abs(outer(0:(nt-1),0:(nt-1),"-"))))))
    dat$y  <-  with(dat, 5 + blockeff + standeff + resid)+rnorm(nrow(dat), 0, .1)
    dat$Time  <-  factor(dat$Time)
    ## saveRDS(dat, file="../../inst/test_data/cov_struct_order.rds",version=2)
}

test_that("cov_struct_order", {
    skip_on_cran()

    fit1  <-  glmmTMB(y ~ (1|Block) + (1|Stand)+ ar1(Time +0|Stand), data = dat)
    expect_equal(unname(fit1$fit$par),
		c(4.98852432, -2.11104196068295, -0.76452645, -0.24762133,  0.08879302,  1.00022657), tol=1e-3)
})

test_that("hom vs het diag", {
    fmhomdiag   <- glmmTMB(Reaction ~ Days + homdiag(Days| Subject), sleepstudy)
    expect_equal(c(VarCorr(fmhomdiag)$cond$Subject),
                 c(69.4182616453357, 0, 0, 69.4182616453357),
                 ## tolerance loosened for var -> SD reparameterization
                 tolerance = 2e-4)

})

test_that("basic hetar1", {
    skip_on_cran()
    vv <- VarCorr(fm_hetar1)[["cond"]]
    expect_equal(attr(vv[[1]], "correlation")[2,1], 0.8861108, tolerance = 1e-5)
    expect_equal(
      attr(vv[[1]], "stddev"),
      c(
        fDays1 = 3.41933, fDays2 = 3.26796, fDays3 = 3.46971, fDays4 = 3.51311, 
        fDays5 = 3.64887, fDays6 = 3.6512, fDays7 = 3.80629, fDays8 = 3.74593,
        fDays9 = 3.80094, fDays10 = 3.54704
      ),
      tolerance = 1e-3
    )
    expect_equal(
      fixef(fm_hetar1)$cond,
      c("(Intercept)" = -0.4591782),
      tolerance = 1e-5
    )
})
