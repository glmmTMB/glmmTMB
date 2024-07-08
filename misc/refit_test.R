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


## example from mebrooks
 
(m1 <- glmmTMB(count ~ mined + (1|site),
   family=poisson, data=Salamanders))

Salamanders2 <- transform(Salamanders, count2=2*count)

refit(m1, Salamanders2$count2)

glmmTMB:::fast_refit(m1, Salamanders2$count2)


## Dispersion model
sim1 <- function(nfac=40, nt=100, facsd=0.1, tsd=0.15, mu=0, residsd=1)
{
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
m0 <- glmmTMB(x ~ sd + (1|t), dispformula=~sd, data=dat)

dat2 <- transform(dat, x2=2*x)
refit(m0, dat2$x2)

debug(fast_refit)
refit(m0, dat2$x2, fast = TRUE)
