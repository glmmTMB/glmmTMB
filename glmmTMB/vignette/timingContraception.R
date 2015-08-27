source("timingFuns.R")
data(Contraception,package="mlmRev")
library("lme4")
library("glmmTMB")
library("glmmADMB")
library("plyr")  ## load before dplyr; for ldply()
library("tidyr")
library("dplyr")

## make sure this is run with optimized build of glmmTMB,
## i.e. "make install" rather than "make quick-install/quick-check"
## (or at least document)
nmax <- 40     ## max replications for glmer/glmmTMB
nmaxADMB <- 2  ## max reps for glmmadmb (much slower)
 
## slow enough that I should consider using something with 
## checkpointing instead ...
tmat <- ldply(seq(nmax),getTimes)  
tmatADMB <- ldply(seq(nmaxADMB),getTimes,which="glmmadmb")

## reshape: wide-to-long, add n values
ff <- function(dd,n=seq(nrow(dd))) {
   mutate(dd,n=n) %>%
    gather(pkg,time,-n)
}
tmatContraception <- rbind(ff(tmat),ff(tmatADMB))

save("tmatContraception",file="contraceptionTimings.rda")

