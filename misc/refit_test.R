devtools::load_all("glmmTMB")

library(glmmTMB)
library(waldo)
fm1 <- glmmTMB(count~mined+(1|spp),
               ziformula=~mined,
               data=Salamanders,
               family=nbinom1)

data("sleepstudy", package = "lme4")
fm2 <- glmmTMB(Reaction ~ Days + (Days | Subject), data = sleepstudy)

fm2B <- glmmTMB(Reaction ~ Days + (Days | Subject), data = sleepstudy, se = FALSE)

## FIXME:: up2date to add control
## single parametric bootstrap step: refit with data simulated from original model
s1 <- simulate(fm1, seed = 101)[[1]]
fm1R <- refit(fm1, s1)
fm1RS <- refit(fm1, s1, fast = TRUE, update_start = FALSE)
## needs work!
all.equal(fm1R, fm1RS)
all.equal(fm1R, fm1RS, tolerance = 1e-4)

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
