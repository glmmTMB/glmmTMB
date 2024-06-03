data("sleepstudy", package = "lme4")
fsleepstudy <- transform(sleepstudy,fDays=cut(Days,c(0,3,6,10),right=FALSE),
                         row=factor(seq(nrow(sleepstudy))))
devtools::load_all("~/R/pkgs/glmmTMB/glmmTMB")
fm_ar1 <- glmmTMB(Reaction ~ 1 +
                      (1|Subject) + ar1(row+0| Subject), fsleepstudy,
                  start = list(theta = c(1,2,0)))
c(VarCorr(fm_ar1)$cond[[1]])

f <- function(ss) {
    fm_ar1 <- glmmTMB(Reaction ~ 1 +
                          (1|Subject) + ar1(row+0| Subject), fsleepstudy,
                      start = list(theta = c(1,ss,0)))
    c(VarCorr(fm_ar1)$cond[[1]])
}
svec <- seq(-1, 3, length = 51)
avec <- sapply(svec, f)
plot(svec, avec, log = "y")

## 7.720272e-05 (bad)

## from head dir:
## git checkout 4e8092357
## Rscript misc/test-ar1var.R
## (requires full recompile each time, sigh)
## git bisect start
## git bisect good    (on 4e8092357, 0.399085)
## git checkout master
## git bisect bad
## ... etc

##  git bisect bad
## Bisecting: 0 revisions left to test after this (roughly 0 steps)
## [5939686f39b169b7d5a81512b6a66025f103d598] update cor calcs notebook [skip ci]
## October 24
