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
prnd <- predict(g0, newdata=nd, allow.new.levels=TRUE)
expect_equal( as.numeric(prnd),
              fixef(g0)$cond[1] + fixef(g0)$cond[2] * nd$Days , tol=1e-10)

context("Catch invalid predictions")

g1 <- glmmTMB(Reaction ~ Days + Subject, sleepstudy)
expect_error( predict(g1, nd) )   ## Error: Requires unknown fixed effect 'Subjectnew'

g2 <- glmmTMB(Reaction ~ us(DaysFac | Subject), sleepstudy)
expect_error( predict(g2, nd) )   ## Error: Not OK due to new parameters

g3 <- glmmTMB(Reaction ~ ar1(DaysFac + 0| Subject), sleepstudy)
expect_warning( predict(g3, nd) ) ## OK: AR1 does not introduce new parameters

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
expect_equal(fitted(tmbm1),fitted(tmbm2))
expect_equal(predict(tmbm1),predict(tmbm2))

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
## type="link"
link_pred <- mypred(~period,dd,fixef(g0_zi)$cond,vcov(g0_zi)$cond)
expect_equal(un(predict(g0_zi,newdata=dd,se.fit=TRUE)),
             link_pred)
## type="conditional"
cond_pred <- mypred(~period,dd,fixef(g0_zi)$cond,vcov(g0_zi)$cond,
                    ff$linkinv,ff$mu.eta)
expect_equal(un(predict(g0_zi,newdata=dd,se.fit=TRUE,type="conditional")),
             cond_pred)
## type="zprob"
zprob_pred <- mypred(~period,dd,fixef(g0_zi)$zi,vcov(g0_zi)$zi,
                    ff$linkinv,ff$mu.eta)
expect_equal(un(predict(g0_zi,newdata=dd,se.fit=TRUE,type="zprob")),
             zprob_pred)
## type="response"  (checking fit only)
expect_equal(unname(predict(g0_zi,newdata=dd,se.fit=TRUE,type="response")$fit),
             cond_pred$fit*(1-zprob_pred$fit))
## type="zlink"
zlink_pred <- mypred(~period,dd,fixef(g0_zi)$zi,vcov(g0_zi)$zi)
expect_equal(un(predict(g0_zi,newdata=dd,se.fit=TRUE,type="zlink")),
             zlink_pred)

context("deprecated zitype parameter")
expect_warning(predict(g0_zi,newdata=dd,zitype="zprob"))
    
    
