## on CRAN: lme4, brms, helper packages
## glmmTMB: need *development version* (>0.1.1) for Salamanders data
## INLA: 
##    install.packages("INLA", repos="https://www.math.ntnu.no/inla/R/stable")
##  (may need to install sp from CRAN separately)
## glmmADMB:
##    devtools::install_github("bbolker/glmmadmb")
## glmmTMB:
##    devtools::install_github("glmmTMB/glmmTMB/glmmTMB")

library(glmmTMB)
library(ggplot2); theme_set(theme_bw())
library(knitr)
library(reshape)
library(plyr)

#Set up cluster
library(parallel)
no_cores = detectCores() - 1
cl = makeCluster(no_cores)
clusterEvalQ(cl, {
  library(glmmTMB)
  library(glmmADMB)
  library(lme4)
  library(brms)
  library(INLA)
  library(mgcv)
  tfun = function(x) unname(system.time(capture.output(x))["elapsed"])
})

#The basic model to fit
mod = glmmTMB(count~spp * mined + (1|site), Salamanders, family="nbinom2")

########################################################
#Benchmark with simulated data sets the same as the original

sims1 = lapply(simulate(mod, nsim=100, seed=111), 
               function(count){ 
                 cbind(count, Salamanders[,c('site', 'mined', 'spp')])
               })

times.tmb = do.call(c, parLapply(cl, sims1, function(x){
  tfun(glmmTMB(count~spp * mined + (1|site), x, family="nbinom2"))}))

times.admb = do.call(c, parLapply(cl, sims1, function(x){
  tfun(glmmadmb(count~spp * mined + (1|site), x, family="nbinom2"))}))

times.lme4 = do.call(c, parLapply(cl, sims1, function(x){
  tfun(glmer.nb(count~spp * mined + (1|site), x))}))

times.brms = do.call(c, parLapply(cl, sims1, function(x){
  tfun(brm(count~spp * mined + (1|site), x, family="negbinomial"))}))

times.inla = do.call(c, parLapply(cl, sims1, function(x){
  tfun(inla(count~spp*mined+f(site, model="iid"), family= "nbinomial", data = x))}))

times.gam = do.call(c, parLapply(cl, sims1, function(x){
  tfun(gam(count ~ spp * mined + s(site, bs = "re"), family = nb, method = "ML", data = x))}))

save(times.tmb, times.admb, times.lme4, times.brms, times.inla, times.gam, file="simfit.RData")

simtimes = data.frame(time=c(times.tmb, times.admb, times.lme4, times.brms, times.inla, times.gam),
                      package=rep(c("glmmTMB", "glmmADMB", "lme4", "brms", "INLA", "mgcv"), each=length(times.tmb)),
                      set=rep(1:length(times.tmb), 6))

ggplot(simtimes)+geom_boxplot(aes(x=package, y=time))+scale_y_log10(breaks=c(1,5,10,50,100,500,1000))+
  ylab("estimation time (seconds)")

ggsave("simfit.png", height=4, width=5)

######################################################## 
#Benchmark with data sets of increasing size

n = nrow(Salamanders)
reps = c(1,2,4,6,8,10)
rows = lapply(reps, function(x) rep(1:n, x))
bigdat = lapply(rows, function(x) Salamanders[x,])

times.tmb = do.call(c, parLapply(cl, bigdat, function(x){
  tfun(glmmTMB(count~spp * mined + (1|site), x, family="nbinom2"))}))

times.admb = do.call(c, parLapply(cl, bigdat, function(x){
  tfun(glmmadmb(count~spp * mined + (1|site), x, family="nbinom2"))}))

times.lme4 = do.call(c, parLapply(cl, bigdat, function(x){
  tfun(glmer.nb(count~spp * mined + (1|site), x))}))

times.brms = do.call(c, parLapply(cl, bigdat, function(x){
  tfun(brm(count~spp * mined + (1|site), x, family="negbinomial"))}))

times.inla = do.call(c, parLapply(cl, bigdat, function(x){
  tfun(inla(count~spp*mined+f(site, model="iid"), family= "nbinomial", data = x, num.threads=1))}))
 
times.gam = do.call(c, parLapply(cl, bigdat, function(x){
  tfun(gam(count ~ spp * mined + s(site, bs = "re"), family = nb, method = "ML", data = x))}))

save(times.tmb, times.admb, times.lme4, times.brms, times.inla, times.gam, file="bigfit.RData")

bigtimes = data.frame(time=c(times.tmb, times.admb, times.lme4, times.brms, times.inla, times.gam),
                      package=rep(c("glmmTMB", "glmmADMB", "lme4", "brms", "INLA", "mgcv"), each=length(times.tmb)),
                      nobs=rep(reps*n, 6))

ggplot(bigtimes, aes(x=nobs, y=time, colour=package))+
  geom_point()+geom_smooth(method="lm", se=FALSE)+
  scale_y_log10(breaks=c(1,5,10,50,100,500,1000))+
  xlab("number of observations in data")+ylab("estimation time (seconds)")+scale_x_log10(breaks=unique(bigtimes$nobs))

ggsave("bigfit.png", height=4, width=5)

########################################################
#Benchmark with data sets with more random effect levels

sims1 = lapply(simulate(mod, nsim=512, seed=111), 
               function(count){ 
                 cbind(count, Salamanders[,c('site', 'mined', 'spp')])
               })
               
n = nrow(Salamanders)
nRE=length(unique(Salamanders$site))

reps = c(1,4,16,64,256,512)
bigdat0 = lapply(reps, function(x) do.call(rbind, sims1[1:x]))
bigdat =  lapply(1:length(reps), function(x)  data.frame(bigdat0[[x]], "grp"=paste0(bigdat0[[x]]$site, rep(1:reps[x], each=n)))) 

times.tmb = do.call(c, parLapply(cl, bigdat, function(x){
  tfun(glmmTMB(count~spp * mined + (1|grp), x, family="nbinom2"))}))

times.inla = do.call(c, parLapply(cl, bigdat, function(x){
  tfun(inla(count~spp*mined+f(grp, model="iid"), family= "nbinomial", data = x))}))

reps = c(1,4,16,64)
bigdat0 = lapply(reps, function(x) do.call(rbind, sims1[1:x]))
bigdat =  lapply(1:length(reps), function(x)  data.frame(bigdat0[[x]], "grp"=paste0(bigdat0[[x]]$site, rep(1:reps[x], each=n)))) 

times.admb = do.call(c, parLapply(cl, bigdat, function(x){
  tfun(glmmadmb(count~spp * mined + (1|grp), x, family="nbinom2", extra.args="-ndi 100000"))}))

times.lme4 = do.call(c, parLapply(cl, bigdat, function(x){
  tfun(glmer.nb(count~spp * mined + (1| grp), x))}))

times.gam = do.call(c, parLapply(cl, bigdat, function(x){
  tfun(gam(count ~ spp * mined + s(grp, bs = "re"), family = nb, method = "ML", data = x))}))

reps = c(1,4,16)
bigdat0 = lapply(reps, function(x) do.call(rbind, sims1[1:x]))
bigdat =  lapply(1:length(reps), function(x)  data.frame(bigdat0[[x]], "grp"=paste0(bigdat0[[x]]$site, rep(1:reps[x], each=n)))) 

times.brms = do.call(c, parLapply(cl, bigdat, function(x){
  tfun(brm(count~spp * mined + (1| grp), x, family="negbinomial"))}))


save(times.tmb, times.inla, times.admb, times.lme4, times.brms, times.gam, file="bigfitRE.RData")

bigtimes = data.frame(time=c(times.tmb, times.inla, times.admb, times.lme4, times.gam, times.brms),
                      package=rep(c("glmmTMB", "INLA", "glmmADMB", "lme4", "mgcv", "brms"),
                                  times=c(length(times.tmb), length(times.inla), length(times.admb), length(times.lme4), length(times.gam), length(times.brms))),
                      nRE=c(rep(c(1,4,16,64,256,512)*nRE, 2),
                            rep(c(1,4,16,64)*nRE, 3),
                            rep(c(1,4,16)*nRE,1)))

ggplot(bigtimes, aes(x=nRE, y=time, colour=package))+
  geom_point()+geom_smooth(method="lm", se=FALSE)+
  scale_y_log10(breaks=c(1,5,10,50,100,500,1000,5000,10000))+
  xlab("number of random effect levels (sites)")+ylab("estimation time (seconds)")+scale_x_log10(breaks=unique(bigtimes$nRE))

ggsave("bigfitRE.png", height=4, width=5.5)
