set.seed(101)
n <- 2000
ngrp <- 20
dd <- data.frame(x=rnorm(n),y=rnorm(n),z=rnorm(n),
                 g1=factor(sample(ngrp,size=n,replace=TRUE)),
                 g2=factor(sample(ngrp,size=n,replace=TRUE)))

library(lme4)
library(glmmTMB)
dd$w <- simulate(~1 + (x+y+z|g1) + (x+y+z|g2),
                 newdata=dd,
                 family=gaussian,
                 newparams=list(beta=1,
                                theta=rep(1,2*10),sigma=1))[[1]]

m1 <- lmer(w~1+(x+y+z|g1), data=dd)
VarCorr(m1)
m2 <- glmmTMB(w~1+(x+y+z|g1), data=dd)
VarCorr(m2)$cond

m3 <- glmmTMB(w~1+rr(x+y+z|g1,3), data=dd)
zapsmall(eigen(VarCorr(m3)$cond$g)$values)
m4 <- glmmTMB(w~1+rr(x+y+z|g1,2), data=dd)
zapsmall(eigen(VarCorr(m4)$cond$g)$values)

debug(glmmTMB:::getXReTrms)
m4 <- glmmTMB(w~1+rr(x+y+z|g1,2), data=dd)
m5 <- glmmTMB(w~1+rr(x+y+z|g1,2)+rr(x+y+z|g2,3), data=dd)
glmmTMB(w~1+rr(x+y+z|g1,2)+(x+y+z|g2), data=dd)
glmmTMB(w~1+rr(x+y+z|g1,n=2)+(x+y+z|g2), data=dd)
