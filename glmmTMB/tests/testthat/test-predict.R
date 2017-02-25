stopifnot(require("testthat"),
          require("glmmTMB"))

data(sleepstudy,
     package = "lme4")
sleepstudy <- transform(sleepstudy, DaysFac = factor(cut(Days,2)) )

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
