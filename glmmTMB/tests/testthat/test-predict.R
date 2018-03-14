stopifnot(require("testthat"),
          require("glmmTMB"))

data(sleepstudy,
     package = "lme4")
sleepstudy <- transform(sleepstudy, DaysFac = factor(cut(Days,2)) )
ssNA <- transform(sleepstudy, Days = replace(Days,c(1,27,93,145), NA))
ssNA2 <- transform(sleepstudy, Days = replace(Days,c(2,49), NA))

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
expect_equal(predict(fm),
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




