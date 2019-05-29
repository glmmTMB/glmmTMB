stopifnot(require("testthat"),
          require("glmmTMB"))

data(sleepstudy, package = "lme4")
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

context("Predicting new levels")

g0 <- glmmTMB(Reaction ~ Days + (Days|Subject), sleepstudy)

test_that("manual prediction of pop level pred", {
    prnd <- predict(g0, newdata=nd, allow.new.levels=TRUE)
    expect_equal( as.numeric(prnd),
                 fixef(g0)$cond[1] + fixef(g0)$cond[2] * nd$Days , tol=1e-10)
})

context("Catch invalid predictions")

test_that("new levels of fixed effect factor", {
    g1 <- glmmTMB(Reaction ~ Days + Subject, sleepstudy)
    expect_error( predict(g1, nd),
                 "Prediction is not possible for unknown fixed effects")
})

test_that("new levels in RE term", {
    g2 <- glmmTMB(Reaction ~ us(DaysFac | Subject), sleepstudy)
    expect_error( predict(g2, nd),
                 "Prediction is not possible for terms")
})

test_that("new levels in AR1 (OK)", {
    g3 <- glmmTMB(Reaction ~ ar1(DaysFac + 0| Subject), sleepstudy)
    expect_warning( predict(g3, nd),
                   ## OK: AR1 does not introduce new parameters
                   "Predicting new random effect levels")
})

context("Predict two-column response case")

fm <- glmmTMB( cbind(count,4) ~ mined, family=betabinomial, data=Salamanders)
expect_equal(predict(fm, type="response"),
             c(0.05469247, 0.29269818)[Salamanders$mined] )

context("Prediction with dispformula=~0")
y <- 1:10
f <- glmmTMB(y ~ 1, dispformula=~0)
expect_equal(predict(f), rep(5.5, 10))

context("Handling NA values in predictions")
ss <- sleepstudy

g0_ex <- update(g0, data=ssNA, na.action=na.exclude)
g0_om <- update(g0, data=ssNA, na.action=na.omit)
pp_ex <- predict(g0_ex)
pp_om <- predict(g0_om)
expect_equal(length(pp_ex),nrow(ssNA))
expect_true(all(is.na(pp_ex)==is.na(ssNA$Days)))
expect_equal(length(pp_om),length(na.omit(ssNA$Days)))
expect_true(!any(is.na(pp_om)))

## na.pass
pp_ndNA <- predict(g0,newdata=ssNA)
expect(all(is.na(ssNA$Days)==is.na(pp_ndNA)))
pp_ndNA2 <- predict(g0,newdata=ssNA2)
expect(all(is.na(ssNA2$Days)==is.na(pp_ndNA2)))

## na.omit
pp_ndNA_om <- predict(g0,newdata=ssNA,na.action=na.omit)
expect_equal(length(pp_ndNA_om),sum(complete.cases(ssNA)))

context("prediction with different binomial specs")

tmbm1 <- glmmTMB(cbind(incidence, size - incidence) ~ period + (1 | herd),
                 data = cbpp, family = binomial)
tmbm2 <- update(tmbm1,incidence/size ~ . , weights = size)

test_that("fitted & predicted agree", {
    expect_equal(fitted(tmbm1),fitted(tmbm2))
    expect_equal(predict(tmbm1),predict(tmbm2))
})

context("zero-inflation prediction")

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

context("deprecated zitype parameter")
expect_warning(predict(g0_zi,newdata=dd,zitype="zprob"))
    
    
context("complex bases")
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
    g3 <- glmmTMB(Reaction~(1|Subject) + scale(Days), sleepstudy)
    expect_equal(predict(g3, newdata=nd, allow.new.levels=TRUE),
                 173.83923026, tolerance=1e-5)
})

test_that("complex bases in dispformula", {
    g4A <- glmmTMB(Reaction~1, sleepstudy)
    g4B <- glmmTMB(Reaction~1,
                   disp=~poly(Days,2), sleepstudy)
    expect_equal(predict(g4A, newdata=nd, se.fit=TRUE),
                 list(fit = 298.507945749154, se.fit = 4.18682101029576),
                 tolerance=1e-5)
    expect_equal(predict(g4B, newdata=nd, se.fit=TRUE),
                 list(fit = 283.656705454758, se.fit = 4.74204256781178))
})
