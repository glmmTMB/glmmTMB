devtools::load_all("glmmTMB")

library(glmmTMB)
library(waldo)
fm1 <- glmmTMB(count~mined+(1|spp),
               ziformula=~mined,
               data=Salamanders,
               family=nbinom1)
## single parametric bootstrap step: refit with data simulated from original model
s1 <- simulate(fm1, seed = 101)[[1]]
fm1R <- refit(fm1, s1)
fm1RS <- refit(fm1, s1, fast = TRUE, update_start = FALSE)
## needs work!
all.equal(fm1R, fm1RS)

## why are $.Phi values different??
waldo::compare(fm1R$modelInfo, fm1RS$modelInfo, ignore_formula_env = TRUE)
## why does $fitted get assigned differently?
fm1R$fitted
fm1RS$fitted

nsim <- 100
options(glmmTMB.fast_refit = FALSE)
t1 <- system.time(b1 <- lme4::bootMer(fm1, FUN=function(x) fixef(x)$zi, nsim=nsim, .progress="txt"))

options(glmmTMB.fast_refit = TRUE)
t2 <- system.time(b2 <- lme4::bootMer(fm1, FUN=function(x) fixef(x)$zi, nsim=nsim, .progress="txt"))
t1
t2

if (requireNamespace("boot")) {
    boot.ci(b1,type="perc")
}
