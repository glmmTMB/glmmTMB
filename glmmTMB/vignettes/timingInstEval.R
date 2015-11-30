source("timingFuns.R")
data(InstEval,package="lme4")
library("lme4")
library("glmmTMB")
library("plyr")  ## load before dplyr; for ldply()
library("tidyr")
library("dplyr")

## make sure this is run with optimized build of glmmTMB,
## i.e. "make install" rather than "make quick-install/quick-check"
## (or at least document)
nvals <- seq(0.1,1,by=0.1)
 
form <- y~service+lectage+studage+(1|d)+(1|s)+(1|dept)
tmat <- ldply(nvals,getTimes,basedata=InstEval,
              form=form,family=NULL,which=c("glmmTMB","lmer"),
              .progress="text")

## reshape: wide-to-long, add n values
ff <- function(dd,n=seq(nrow(dd))) {
   mutate(dd,n=nvals) %>%
    gather(pkg,time,-n)
}
tmatInstEval <- ff(tmat)
save("tmatInstEval",file="InstEvalTimings.rda")

