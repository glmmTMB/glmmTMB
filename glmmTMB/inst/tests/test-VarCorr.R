stopifnot(require("testthat"),
          require("glmmTMB"))

context("VarCorr Testing")
##       ---------------

data("Orthodont", package="nlme")
fm1 <- glmmTMB(distance ~ age + (age|Subject), data = Orthodont)

data("Pixel", package="nlme")
nPix <- nrow(Pixel)
if(FALSE) ## segmentation fault !!
fmPix <- lmer(pixel ~ day + I(day^2) + (day | Dog) + (1 | Side/Dog), data = Pixel)
## TMB has received an error from Eigen. The following condition was not met:
## index >= 0 && index < size()
## Please check your matrix-vector bounds etc., or run your program through a debugger.
## " Abgebrochen  (Speicherabzug geschrieben) " -- (core dump)

set.seed(12345)
dd <- data.frame(a=gl(10,100), b=rnorm(1000))
library(lme4)
test2 <- simulate(~1+(b|a), newdata=dd, newparams=list(beta=c(1), theta=c(1,1,1)), family=poisson)
mydata <- cbind(dd, test2)
mydata$sim_1[sample(c(FALSE,TRUE), 1000, prob=c(.3,.7), replace=TRUE)] <- 0
# not simulated this way, but returns right structure
gm <- glmmTMB(sim_1~1+(b|a), zi=~1+(b|a), data=mydata,family=poisson())
gm$obj$env$report()



quit()
##===  Not yet :

(vc <- VarCorr(fm1))  ## default print method: standard dev and corr
## both variance and std.dev.
print(vc,comp=c("Variance","Std.Dev."),digits=2)
## variance only
print(vc,comp=c("Variance"))
as.data.frame(vc)
as.data.frame(vc,order="lower.tri")
