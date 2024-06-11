## this is the test/imports that are weird:

## with current master (HEAD), fresh install:
library(glmmTMB)
sf <- system.file("test_data/models.rda", package = "glmmTMB")
L <- load(sf)
try(fm_ar1) ## incorrect number of args
get_theta <- function(m) {
    ee <- m$obj$env
    pp <- ee$last.par.best
    split(pp, names(pp))$theta
}
get_theta(fm_ar1) ## OK (-0.4593) ## now -4.772030 ... ????
fm_ar2 <- up2date(fm_ar1)
get_theta(fm_ar2)

rm(fm_ar1)
file_ok <- gt_load("test_data/models.rda")
get_theta(fm_ar1) ## OK

data("sleepstudy", package = "lme4")
fsleepstudy <- transform(sleepstudy,fDays=cut(Days,c(0,3,6,10),right=FALSE),
                         row=factor(seq(nrow(sleepstudy))))

library(glmmTMB)
fm_ar1 <- glmmTMB(Reaction ~ 1 +
                      (1|Subject) + ar1(row+0| Subject), fsleepstudy)
ee <- fm_ar1$obj$env
with(ee, last.par.best[-random])

pars1 <- c(301.3836838,4.3687610,-0.4593243,4.1074652,1.7899137)
pars2 <- c(301.3836838,2.184327 ,-4.862984 ,4.1074652,1.7899137)
fm_ar1$obj$fn(pars1)
fm_ar1$obj$fn(pars2) ## ??
##       beta        betad      theta      theta      theta 
## 301.382851    2.184327  -4.862984   4.107471   1.789900
##        beta       betad       theta       theta       theta 
## 301.3836838   4.3687610  -0.4593243   4.1074652   1.7899137
## almost identical LLs (?) with values of `betad`, `theta[1]` shifted

## 'log Lik.' -890.1838 (df=5)
## 'log Lik.' -890.1844 (df=5)

## check 2D grid?

## re-install with older (CppAD)

## git checkout reparam_gauss_family; re-install

## ar1 model was re-run here:
## commit 9237e5eff01babb788ba5aed8ecd3607398922ff
## Author: Ben Bolker <bbolker@gmail.com>
## Date:   Mon Oct 30 11:39:42 2023 -0400
##     first pass at Gaussian var -> SD param shift


## with SD param:


library(glmmTMB)
sf <- system.file("test_data/models.rda", package = "glmmTMB")
L <- load(sf)
sigma(fm2)
## should be modifying par, parfull
debug(glmmTMB::up2date)
fm2B <- up2date(fm2, update_gauss_disp = TRUE)
sigma(fm2B)
glmmTMB:::getParList(fm2)$betad
