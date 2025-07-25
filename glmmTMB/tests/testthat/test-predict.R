stopifnot(require("testthat"),
          require("glmmTMB"))

sleepstudy <- transform(sleepstudy, DaysFac = factor(cut(Days,2)) )
ssNA <- transform(sleepstudy, Days = replace(Days,c(1,27,93,145), NA))
ssNA2 <- transform(sleepstudy, Days = replace(Days,c(2,49), NA))

data(cbpp, package = "lme4")
set.seed(101)
cbpp_zi <- cbpp
cbpp_zi[sample(nrow(cbpp),size=15,replace=FALSE),"incidence"] <- 0

## 'newdata'
nd <- subset(sleepstudy, Subject=="308", select=-1)
nd$Subject <- "new"
nd$DaysFac <- "new"

test_that("manual prediction of pop level pred", {
    prnd <- predict(fm2, newdata=nd, allow.new.levels=TRUE)
    expect_equal( as.numeric(prnd),
                 fixef(fm2)$cond[1] + fixef(fm2)$cond[2] * nd$Days , tol=1e-10)
})

test_that("population-level prediction", {
    prnd <- predict(fm2)
    expect_equal(length(unique(prnd)),180)
    prnd2 <- predict(fm2, re.form=~0)
    prnd3 <- predict(fm2, re.form=NA)
    expect_equal(prnd2,prnd3)
    expect_equal(length(unique(prnd2)),10)
    ## make sure we haven't messed up any internal structures ...
    prnd4 <- predict(fm2)
    expect_equal(prnd, prnd4)
})

test_that("new levels of fixed effect factor", {
    skip_on_cran()
    g1 <- glmmTMB(Reaction ~ Days + Subject, sleepstudy)
    expect_error( predict(g1, nd),
                 "Prediction is not possible for unknown fixed effects")
})

test_that("new levels in RE term", {
    skip_on_cran()
    suppressWarnings(g2 <- glmmTMB(Reaction ~ us(DaysFac | Subject), sleepstudy))
    expect_error( predict(g2, nd),
                 "Prediction is not possible for terms")
})

test_that("new levels in AR1 (OK)", {
    skip_on_cran()
    g3 <- glmmTMB(Reaction ~ ar1(DaysFac + 0| Subject), sleepstudy)
    expect_warning( predict(g3, nd),
                   ## OK: AR1 does not introduce new parameters
                   "Predicting new random effect levels")
})

test_that("two-column response", {
    skip_on_cran()
    fm <- glmmTMB( cbind(count,4) ~ mined, family=betabinomial,
                  data=Salamanders)
    expect_equal(predict(fm, type="response"),
                 c(0.05469247, 0.29269818)[Salamanders$mined] )
})

test_that("Prediction with dispformula=~0", {
    skip_on_cran()
    y <- 1:10
    f <- glmmTMB(y ~ 1, dispformula=~0, data = NULL)
    expect_equal(predict(f), rep(5.5, 10))
})

ss <- sleepstudy

fm2_ex <- update(fm2, data=ssNA, na.action=na.exclude)
fm2_om <- update(fm2, data=ssNA, na.action=na.omit)
pp_ex <- predict(fm2_ex)
pp_om <- predict(fm2_om)

test_that("NA values in predictions", {
    expect_equal(length(pp_ex),nrow(ssNA))
    expect_true(all(is.na(pp_ex)==is.na(ssNA$Days)))
    expect_equal(length(pp_om),length(na.omit(ssNA$Days)))
    expect_true(!any(is.na(pp_om)))
})

## na.pass
test_that("na.pass", {
    pp_ndNA <- predict(fm2,newdata=ssNA)
    expect(all(is.na(ssNA$Days)==is.na(pp_ndNA)),
           failure_message="NAs don't match with na.pass+predict")
    pp_ndNA2 <- predict(fm2,newdata=ssNA2)
    expect(all(is.na(ssNA2$Days)==is.na(pp_ndNA2)),
           failure_message="NAs don't match with na.pass+predict+newdata")
})

## na.omit
test_that("na.omit", {
    pp_ndNA_om <- predict(fm2,newdata=ssNA,na.action=na.omit)
    expect_equal(length(pp_ndNA_om),sum(complete.cases(ssNA)))
})

tmbm1 <- glmmTMB(cbind(incidence, size - incidence) ~ period + (1 | herd),
                 data = cbpp, family = binomial)
tmbm2 <- update(tmbm1,incidence/size ~ . , weights = size)

test_that("different binomial specs: fitted & predicted agree", {
    expect_equal(fitted(tmbm1),fitted(tmbm2))
    expect_equal(predict(tmbm1),predict(tmbm2))
})

## context("zero-inflation prediction")

g0_zi <- update(tmbm2, ziformula = ~period)
un <- function(x) lapply(x,unname)
mypred <- function(form,dd,cc,vv,linkinv=identity,mu.eta=NULL) {
    X <- model.matrix(form,dd)
    pred <- drop(X %*% cc)
    se <- drop(sqrt(diag(X %*% vv %*% t(X))))
    if (!is.null(mu.eta)) se <- se*mu.eta(pred)
    pred <- linkinv(pred)
    return(un(list(fit=pred,se.fit=se)))
}
## FIXME: predictions should have row names of data
dd <- data.frame(unique(cbpp["period"]),size=1,herd=NA)
ff <- make.link("logit")

test_that("type='link'", {
    link_pred <- mypred(~period,dd,fixef(g0_zi)$cond,vcov(g0_zi)$cond)
    expect_equal(un(predict(g0_zi,newdata=dd,se.fit=TRUE)),
                 link_pred)
})

test_that("various types", {
    cond_pred <- mypred(~period,dd,fixef(g0_zi)$cond,vcov(g0_zi)$cond,
                        ff$linkinv,ff$mu.eta)
    expect_equal(un(predict(g0_zi,newdata=dd,se.fit=TRUE,type="conditional")),
                 cond_pred)
    zprob_pred <- mypred(~period,dd,fixef(g0_zi)$zi,vcov(g0_zi)$zi,
                         ff$linkinv,ff$mu.eta)
    expect_equal(un(predict(g0_zi,newdata=dd,se.fit=TRUE,type="zprob")),
                 zprob_pred)
    expect_equal(unname(predict(g0_zi,newdata=dd,se.fit=TRUE,type="response")$fit),
                 cond_pred$fit*(1-zprob_pred$fit))
})

test_that("type='zlink'", {
    zlink_pred <- mypred(~period,dd,fixef(g0_zi)$zi,vcov(g0_zi)$zi)
    expect_equal(un(predict(g0_zi,newdata=dd,se.fit=TRUE,type="zlink")),
                 zlink_pred)
})

test_that("deprecated zitype parameter", {
    expect_warning(predict(g0_zi,newdata=dd,zitype="zprob"))
})


## context("complex bases")
data("sleepstudy",package="lme4")
nd <- data.frame(Days=0,
                 Subject=factor("309", levels=levels(sleepstudy$Subject)))

test_that("poly", {
    g1 <- glmmTMB(Reaction~poly(Days,3), sleepstudy)
    expect_equal(predict(g1, newdata=data.frame(Days=0)),
                 255.7690496, tolerance=1e-5)
})
test_that("splines", {
    if (getRversion()>="3.5.1") {
        ## work around predict/predvars bug in 3.5.0 & previous versions
        g2 <- glmmTMB(Reaction~splines::ns(Days,5), sleepstudy)
    } else {
        library(splines)
        g2 <- glmmTMB(Reaction~ns(Days,5), sleepstudy)
    }
    expect_equal(predict(g2, newdata=data.frame(Days=0)),257.42672,
                 tolerance=1e-5)
})

test_that("scale", {
    skip_on_cran()
    g3 <- glmmTMB(Reaction~scale(Days), sleepstudy)
    expect_equal(predict(g3, newdata=data.frame(Days=0)),
                 251.40507651, tolerance=1e-5)
})
test_that("poly_RE", {
    g1 <- glmmTMB(Reaction~(1|Subject) + poly(Days,3), sleepstudy)
    expect_equal(predict(g1, newdata=nd, allow.new.levels=TRUE),
                 178.1629812, tolerance=1e-5)
})
test_that("splines_RE", {
    if (getRversion()>="3.5.1") {
        g2 <- glmmTMB(Reaction~(1|Subject) + splines::ns(Days,5), sleepstudy)
    } else {
        library(splines)
        g2 <- glmmTMB(Reaction~(1|Subject) + ns(Days,5), sleepstudy)
    }
    expect_equal(predict(g2, newdata=nd, allow.new.levels=TRUE),
                 179.7784754, tolerance=1e-5)
})

test_that("scale_RE", {
    skip_on_cran()
    g3 <- glmmTMB(Reaction~(1|Subject) + scale(Days), sleepstudy)
    expect_equal(predict(g3, newdata=nd, allow.new.levels=TRUE),
                 173.83923026, tolerance=1e-5)
})

test_that("complex bases in dispformula", {
    skip_on_cran()
    g4A <- glmmTMB(Reaction~1, sleepstudy)
    g4B <- glmmTMB(Reaction~1,
                   disp=~poly(Days,2), sleepstudy)
    expect_equal(predict(g4A, newdata=nd, se.fit=TRUE),
                 list(fit = 298.507945749154, se.fit = 4.18682101029576),
                 tolerance=1e-5)
    expect_equal(predict(g4B, newdata=nd, se.fit=TRUE),
                 list(fit = 283.656705454758, se.fit = 4.74204256781178),
                 tolerance = 1e-6)
})

test_that("fix_predvars works for I(x^2)", {
    skip_on_cran()
    ## GH512; @strengejacke
    set.seed(123)
    n <- 500
    d <- data.frame(
        y = rbinom(n, size = 1, prob = .2),
        x = rnorm(n),
        site = sample(letters, size = n, replace = TRUE),
        area = sample(LETTERS[1:9], size = n, replace = TRUE)
    )
    form <- y ~ x + I(x^2) + I(x^3) + (1 | area)
    m1 <- lme4::glmer(form, family = binomial("logit"), data = d)
    m2 <- glmmTMB(form, family = binomial("logit"), data = d)
    nd <- data.frame(x = c(-2, -1, 0, 1, 2), area = NA)
    p1 <- predict(m1, newdata = nd, type = "link", re.form = NA)
    p2 <- predict(m2, newdata = nd, type = "link")
    expect_equal(unname(p1),unname(p2), tolerance=1e-4)
})

test_that("more predvars stuff (I()) (GH #853)", {
    set.seed(100)
    N <- 100
    x1 <- rnorm(N)
    Y <- rpois(N, lambda = exp(x1))
    df <- data.frame(Y=Y, x1=x1)

    ##Base model
    mod <- suppressWarnings(glmmTMB(Y ~ x1 + I((x1+10)^2) + I((x1+10)^3),
                                    data = df, family = "poisson"))
    p1 <- predict(mod, newdata = df, type = "response")
    expect_equal(fitted(mod), predict(mod, newdata = df, type = "response"))
})

test_that("predvars with different ns() in fixed and disp (GH #845)", {
    library(splines)
    x <- glmmTMB(
        am ~ ns(wt, df = 3), 
        dispformula = ~ ns(wt, df = 2), 
        data = mtcars
    )
    newdata <- data.frame(
        wt = seq(min(mtcars$wt), max(mtcars$wt), length.out = 3)
    )
    expect_equal(predict(x, newdata = newdata),
                 c(1.00149139390868, 0.367732526652086, 9.21516947505197e-06))
})

test_that("predvars with differing splines in fixed and RE (GH#632)", {
    library(splines)
    data(sleepstudy,package="lme4")
    form <- Reaction ~ ns(Days, df = 3) + (ns(Days, df = 2)|Subject)
    ## need better/non-default starting values after var -> SD param change
    m4 <- glmmTMB(form, data = sleepstudy,
                  start = list(theta = c(3.28, 4.34, 3.79, 0, 0, 0)))
    m4B <- lme4::lmer(form, data = sleepstudy, REML = FALSE)
    nd <- data.frame(Days = 4:6, Subject = "372")
    predict(m4B, newdata = nd, re.form = NULL, type = "response")
    pp <- predict(m4, newdata = nd, re.form = NULL, type = "response")
    ppB <- predict(m4B, newdata = nd, re.form = NULL, type = "response")
    ## plot(Reaction ~ Days, data = subset(sleepstudy, Subject == "372"))
    ## points(4:6, pp, col = 2, pch = 16)
    ## ?? results change on var to SD switch for Gaussian model?
    expect_equal(pp, c(309.103652912868, 321.193466901353, 333.568337949647),
                 tolerance = 1e-4)
    expect_equal(unname(ppB), pp, tolerance = 1e-4)
})

test_that("contrasts carried over", {
    skip_on_cran()
    ## GH 439, @cvoeten
    iris2 <- transform(iris,
                       grp=factor(c("a","b")))
    contrasts(iris2$Species) <- contr.sum
    contrasts(iris2$grp) <- contr.sum
    mod1 <- glmmTMB(Sepal.Length ~ Species, iris)
    mod2 <- glmmTMB(Sepal.Length ~ Species, iris2)
    ## create prediction frame with a new species
    iris3 <- iris[1,]
    iris3$Species <- "extra"
    ## these are not *exactly* equal because of numeric differences
    ##  when estimating parameters differently ... (?)
    expect_equal(predict(mod1), predict(mod2), tolerance=1e-6)
    ## make sure we actually imposed contrasts correctly/differently
    expect_false(isTRUE(all.equal(fixef(mod1)$cond,fixef(mod2)$cond)))
    expect_warning(predict(mod1,newdata=iris2), "contrasts mismatch")
    expect_equal(predict(mod1,newdata=iris2,allow.new.levels=TRUE),
                 predict(mod1,newdata=iris))
    mod3 <- glmmTMB(Sepal.Length ~ 1|Species, iris)
    expect_equal(c(predict(mod3,newdata=data.frame(Species="ABC"),
                           allow.new.levels=TRUE)),
                 5.843333, tolerance=1e-6)
    mod4 <- glmmTMB(Sepal.Length ~ grp + (1|Species), iris2)
    expect_equal(c(predict(mod4, newdata=data.frame(Species="ABC",grp="a"),
                           allow.new.levels=TRUE)),
                 5.839998, tolerance=1e-6)
    ## works with char rather than factor in new group vble
    expect_equal(predict(mod3, newdata=iris3, allow.new.levels=TRUE),
                 5.843333, tolerance=1e-6)
    zipm3 <- glmmTMB(
        count ~ spp * mined + (1 | site:mined),
        Salamanders,
        family = "poisson"
    )
    ss <- transform(Salamanders, site = factor(site, ordered = FALSE))
    zz <- update(zipm3, data = ss)
    expect_no_warning(p1 <- predict(zipm3, newdata = head(Salamanders, 22)))
    p2 <- predict(zz, newdata = head(Salamanders, 22))
    p3 <- predict(zz, newdata = head(ss, 22))
    p4 <- predict(zipm3, newdata = head(Salamanders, 23))
    expect_equal(p1, p2)
    expect_equal(p2, p3)
    expect_equal(p3, p4[1:22])
})

test_that("dispersion", {
    mod5 <- glmmTMB(Sepal.Length ~ Species, disp=~ Species, iris)
    expect_equal(length(unique(predict(mod5, type="disp"))), length(unique(iris$Species)))
    expect_equal(length(unique(predict(mod5, type="disp", se.fit=TRUE)$se.fit)), length(unique(iris$Species)))
})

test_that("offset-only model (GH #625)", {
    skip_on_cran()
    owls_nb0 <- glmmTMB(SiblingNegotiation ~ offset(log(BroodSize)),
                        family = nbinom2(),
                        data=Owls)
    expect_equal(mean(predict(owls_nb0)), 1.88220473712677)
})

test_that("fast prediction", {
    ## use tighter-than-default tolerances
    ##
    expect_equal(predict(fm2,fast=FALSE),predict(fm2,fast=TRUE), tolerance=1e-13)
    expect_equal(predict(fm2, type="response",fast=FALSE),
                 predict(fm2, type="response", fast=TRUE),
                 tolerance=1e-13)
    ## handling NAs etc.
    expect_equal(pp_ex, predict(fm2_ex, fast=FALSE))
})

test_that("inverse-link prediction", {
  skip_on_cran()
  ## example from John Maindonald (GH #696)
  ## this highlights a particular case where the prediction on the (cloglog) link scale
  ## is large (3.98), which leads to a prediction of 1.0 unless the cloglog-inverse-link
  ## function is clamped (as in make.link("cloglog")'s version)
  ffly <- read.csv(system.file("test_data", "ffly.csv", package="glmmTMB"))
  ffly$obs <- factor(ffly$obs)
  form1 <- cbind(Dead,Live)~0+trtGp/TrtTime+(1|obs)+(1|trtGpRep)
  ObsTMB.cll <- glmmTMB(form1, family=binomial(link="cloglog"), data=ffly)
  p0 <- predict(ObsTMB.cll, re.form=NA)[63]
  p0R <- make.link("cloglog")$linkinv(p0)
  p1 <- predict(ObsTMB.cll, re.form=NA, type="response")[63]
  expect_equal(p0R, p1)
})

test_that("fast prediction not allowed with NA (correct errors)", {
  expect_error(predict(fm2, re.form=NA, fast=TRUE),
               "fast=TRUE is not compatible")
  expect_equal(predict(fm2, re.form=NA, fast=FALSE),
               predict(fm2, re.form=NA, fast=NULL))
})

test_that("zlink/zprob return appropriate values with non-ZI model (GH#798)", {
  p1 <- predict(fm2, type = "zlink")
  expect_equal(length(p1), nrow(sleepstudy))
  expect_true(all(p1 == -Inf))
  p2 <- predict(fm2, type = "zprob")
  expect_equal(length(p2), nrow(sleepstudy))
  expect_true(all(p2 == 0))
})

test_that("correct conditional/response predictions for truncated distributions", {
    set.seed(42)
    N <- 100
    df <- data.frame(p1 = rpois(N, 1),
                     nb1 = rnbinom(N, mu = 1, size = 1),
                     x = rnorm(N)) |>
        transform(
            ## zero-inflated versions
            zp1 = p1 * rbinom(N, size = 1, prob = 0.5),
            znb1 = nb1 * rbinom(N, size = 1, prob = 0.5),
            ## truncated versions (NAs will be dropped)
            tp1 = ifelse(p1 == 0, NA, p1),
            tnb1 = ifelse(nb1 == 0, NA, nb1))

    f_zp1 <- glmmTMB(zp1 ~ x,
                     zi= ~ 1,
                     family=truncated_poisson(link="log"),
                     data=df)

    f_znb1 <- update(f_zp1, znb1 ~ .,
                     family = truncated_nbinom1)

    f_znb2 <- update(f_zp1, znb1 ~ .,
                     family = truncated_nbinom2)


    testfun <- function(model, response, distrib) {
        zp1 <- predict(model, type="zprob")
        cm1 <- predict(model, type="conditional")
        mu1 <- predict(model, type="response")

        ## compute zero-trunc by hand
        eta <- predict(model, type = "link")
        cm2 <- exp(eta)/(1-distrib(0, exp(eta), sigma(model)))
        expect_equal(cm1, cm2)

        expect_equal(mu1, cm1*(1-zp1))

    }

    ## versions of distrib functions that can be plugged into testfun()
    my_dpois <- function(x, lambda, ...) dpois(x, lambda)
    my_nb2 <- function(x, mu, size) dnbinom(x, mu = mu, size = size)
    my_nb1 <- function(x, mu, phi) {
        ## var = mu*(1+mu/k) = mu*(1+phi) -> phi = mu/k -> k = mu/phi
        dnbinom(x, mu = mu, size = mu/phi)
    }

    testfun(f_zp1, "zp1", my_dpois)
    testfun(f_znb2, "znb1", my_nb2)
    testfun(f_znb1, "znb1", my_nb1)

})

test_that("predict warns about ignored args", {
    expect_warning(predict(fm2, bad_args = TRUE), "bad_args")
})

## GH #873
test_that("nzprob doesn't segfault", {
    skip_on_cran()
    model2 <- glmmTMB(
        count ~ cover + mined + (1 | site),
        ziformula = ~ cover + mined,
        family = truncated_poisson(),
        data = Salamanders
    )
    pp <- stats::predict(
               model2,
               newdata = Salamanders,
               type = "link",
               re.form = NULL,
               allow.new.levels = FALSE
               )
    expect_equal(head(pp, 3),
                 c(0.465946249085321, 0.206712238705304, 0.133580349579438))
})

## GH #873 continued
test_that("nzprob computed for non-fast pred", {
    set.seed(101)
    dd <- data.frame(y = rpois(5, lambda = 1))
    m1 <- glmmTMB(
        y ~ 1,
        ziformula = ~ 1,
        data = dd,
        family = truncated_poisson()
    )
    expect_identical(predict(m1, type = "response"),
                     predict(m1, type = "response", fast = FALSE))
    ## non-pos-def Hessian, ignore
    m2 <- suppressWarnings(update(m1, family = truncated_nbinom1))
    expect_identical(predict(m2, type = "response"),
                     predict(m2, type = "response", fast = FALSE))
    ## non-pos-def Hessian, ignore
    m2 <- suppressWarnings(update(m1, family = truncated_nbinom2))
    ## need more data to fit compois, genpois
    dd2 <- data.frame(y = rpois(100, lambda = 1))
    m2 <- update(m1, family = truncated_compois, data = dd2)
    expect_identical(predict(m2, type = "response"),
                     predict(m2, type = "response", fast = FALSE))
    ## suppress NA/NaN function eval warning
    m2 <- suppressWarnings(update(m1, family = truncated_genpois, data = dd2))
    expect_identical(predict(m2, type = "response"),
                     predict(m2, type = "response", fast = FALSE))
})

test_that("pop-level prediction with missing grouping vars (GH #923)",
{
    fm20 <- glmmTMB(Reaction ~ 1 + (Days|Subject), sleepstudy)
    predict(fm20, re.form = NA, newdata = data.frame(matrix(ncol = 0, nrow=length(sleepstudy))))
    expect_equal(predict(fm2, re.form = NA),
                 predict(fm2, newdata = sleepstudy[c("Days")], re.form = NA))
    ## suppress non-pos-def Hessian warning, this is a silly example
    fmnasty <- suppressWarnings(
        glmmTMB(Reaction ~ 1 + (log(Days+1)|Subject), sleepstudy)
    )
    expect_equal(predict(fmnasty, re.form = NA),
                 predict(fmnasty, newdata = data.frame(matrix(ncol=0, nrow = nrow(sleepstudy))), re.form = NA))
})

    
test_that("weights with attributes are OK", {
    data(iris)
    d <- as.data.frame(expand.grid(
        Species = unique(iris$Species),
        Petal.Width = 2,
        wg = NA
    ))
    
    set.seed(101)
    iris$wg <- abs(rnorm(nrow(iris), 1, 0.1))
    
    attr(iris$wg, "label") <- "weighting variable"
    
    m <- glmmTMB(
        Petal.Length ~ Petal.Width + (1 | Species),
        weights = wg,
        data = iris
    )
    p <- predict(m, newdata = d, re.form = NULL)
    expect_equal(head(p),
                 c(3.35535960506176, 4.98184089632094, 5.50821757779119))
})

test_that("pearson resids of ZI models", {
    if (requireNamespace("pscl")) {
        data(bioChemists, package = "pscl")
        m <- glmmTMB::glmmTMB(
                          art ~ fem + mar + kid5 + ment,
                          ziformula = ~ kid5 + phd,
                          family = poisson(),
                          data = bioChemists
                      )
        expect_equal(sum(residuals(m, type = "pearson")^2),
                     1258.55399982061)

    }
})

## GH 1189
test_that("allow.new.levels TRUE when re.form = NA", 
          {
              nd <- subset(sleepstudy, Subject=="308", select=-1)[1,]
              nd$Subject <- "new"
              
              expect_warning(predict(fm2, newdata=nd),
                           "Predicting new random effect levels")
              expect_equal(predict(fm2, newdata = nd, re.form = NA),
                           251.404341632295,
                           tolerance = 1e-6)

              suppressWarnings(g1 <- glmmTMB(Reaction ~ 1 + (1 | Subject/Days),
                                             sleepstudy))
              expect_equal(predict(g1,
                                   newdata = data.frame(Days = NA, Subject = NA),
                                   re.form = NA), 298.507889474305,
                           tolerance = 1e-6)
})

## https://stackoverflow.com/q/77517125/190277
test_that("prediction from rank-deficient X matrices", {
  x <- c("A", "B", "C", "D"); y <- c("exposed", "ref1", "ref2")
  set.seed(123)
  dat <- data.frame(time = rep(x, each=20, times=3),
                    lake = rep(y, each = 80),
                    min = runif(n=240, min=4.5, max=5.5),
                    count = rnbinom(n=240,mu=10,size=100))
  dat2 <- subset(dat, time!="A" | lake !="ref1")
  suppressMessages(
    model <-glmmTMB(count~time*lake, family=nbinom1,
                    control = glmmTMBControl(rank_check = "adjust"),
                    offset=log(min), data=dat2))
  pred_data <- data.frame(
    lake = rep(c("exposed", "ref1", "ref2"), c(4L, 3L, 4L)),
    min = 1,
    time = c(x, x[-1], x)  ## drop level A for lake 2
  )

  pp <- suppressMessages(
    predict(model, newdata = pred_data, type = 'response', se.fit=FALSE)
  )
  expect_equal(head(pp, 3), c(1.89088787010282, 1.83544974589453, 2.02668291828295),
               tolerance = 1e-6)
})
             

