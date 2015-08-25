stopifnot(require("testthat"),
          require("glmmTMB"),
          require("lme4"))

context("VarCorr Testing")
##       ---------------

data("Orthodont", package="nlme")
fm1 <- glmmTMB(distance ~ age + (age|Subject), data = Orthodont)
lm1 <-    lmer(distance ~ age + (age|Subject), data = Orthodont) # to compare
fm1.r <- fm1$obj$env$report()
fm1.r$corr[[1]] ## very highly correlated - ?

data("Pixel", package="nlme")
## nPix <- nrow(Pixel)
if(FALSE) ## segmentation fault !!
fmPix1 <- glmmTMB(pixel ~ day + I(day^2) + (day | Dog) + (1 | Dog/Side), data = Pixel)
## no segfault; just 13 warnings:
unique(warnings()) ## 3 different ones

fmPix2 <- glmmTMB(pixel ~ day + I(day^2) + (day | Dog) + (1 | Side/Dog), data = Pixel)
## Warning messages:
## 1: In mapply(parFun, covCode, blksize, SIMPLIFY = FALSE) :
##   longer argument not a multiple of length of shorter
## 2: In mapply(list, blockReps = nreps, blockSize = blksize, blockNumTheta = blockNumTheta,  :
##   longer argument not a multiple of length of shorter

## "manual"  (1 | Dog / Side) :
fmPix3 <- glmmTMB(pixel ~ day + I(day^2) + (day | Dog) + (1 | Dog) + (1 | Side:Dog), data = Pixel)
## 11 warnings
unique(warnings()) ## all 11 are :
## In nlminb(start = par, objective = fn, gradient = gr) : NA/NaN function evaluation

## "manual"  (1 | Side / Dog) :
fmPix4 <- glmmTMB(pixel ~ day + I(day^2) + (day | Dog) + (1 | Side) + (1 | Side:Dog), data = Pixel)
## no warnings... good!
str(fmP4.r <- fmPix4$obj$env$report())
## List of 4
##  $ corrzi: list()
##  $ sdzi  : list()
##  $ corr  :List of 3
##   ..$ : num [1, 1] 1
##   ..$ : num [1:2, 1:2] 1 -1 -1 1
##   ..$ : num [1, 1] 1
##  $ sd    :List of 3
##   ..$ : num 28.9
##   ..$ : num [1:2] 1.53 1.51
##   ..$ : num 0.000165
fmP4.r $ corr

set.seed(12345)
dd <- data.frame(a=gl(10,100), b = rnorm(1000))
test2 <- simulate(~1+(b|a), newdata=dd, family=poisson,
                  newparams= list(beta = c("(Intercept)" = 1),
                                  theta = c(1,1,1)))
## Zero-inflation : set all i.0 indices to 0:
i.0 <- sample(c(FALSE,TRUE), 1000, prob=c(.3,.7), replace=TRUE)
test2[i.0, 1] <- 0
str(mydata <- cbind(dd, test2))
## The zeros in the 10 groups:
xtabs(~ a + (sim_1 == 0), mydata)


# not simulated this way, but returns right structure
gm <- glmmTMB(sim_1 ~ 1+(b|a), zi = ~1+(b|a), data=mydata, family=poisson())
## eight updateCholesky() warnings .. which will suppress *unless* they are in the last iter.
str(gm.r <- gm$obj$env$report())
## List of 4
##  $ corrzi:List of 1
##   ..$ : num [1:2, 1:2] 1 0.929 0.929 1
##  $ sdzi  :List of 1
##   ..$ : num [1:2] 3.03e-05 1.87e-04
##  $ corr  :List of 1
##   ..$ : num [1:2, 1:2] 1 0.921 0.921 1
##  $ sd    :List of 1
##   ..$ : num [1:2] 0.779 1.575



quit()
##===  Not yet :

(vc <- VarCorr(fm1))  ## default print method: standard dev and corr
## both variance and std.dev.
print(vc,comp=c("Variance","Std.Dev."),digits=2)
## variance only
print(vc,comp=c("Variance"))
as.data.frame(vc)
as.data.frame(vc,order="lower.tri")
