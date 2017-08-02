## Demonstrate optimization issue using cloglog example from
## https://github.com/bbolker/mixedmodels-misc/blob/master/notes/cloglogsim.rmd

library(Matrix)
library(glmmTMB)
library(MASS)
artSim <- function(){
    ##
    ## Function to simulate "artificial" data which is at least superficially
    ## similar to some real data.
    ##
    link    <- "cloglog"
    B       <- binomial(link=link)
    linkfun <- B$linkfun
    linkinv <- B$linkinv

    ## Construct (artificial) treatment factor, covariate, and
    ## (random) replicate factor.
    x    <- seq(0,28,by=2)
    Trt  <- LETTERS[1:24]
    Rep  <- 1:3 ## Three reps per treatment.
    Xdat <- expand.grid(x=x,Trt=Trt,Rep=Rep)
    uRep <- with(Xdat,factor(paste0(Rep,Trt)))
    Xdat$Rep <- with(Xdat,factor(as.numeric(uRep)))

    beta0 <- seq(-3,0.45,by=0.15)
    beta1 <- rep(seq(0.05,0.3,by=0.05),4)
    names(beta0) <- Trt
    names(beta1) <- Trt
    Sigma <- matrix(c(0.06,-0.001,-0.001,0.0001),nrow=2)

    lb0   <- beta0[match(Xdat$Trt,names(beta0))]
    lb1   <- beta1[match(Xdat$Trt,names(beta1))]
    nrep  <- 72
    imat  <- match(Xdat$Rep,1:nrep)
    Z     <- mvrnorm(nrep,c(0,0),Sigma)[imat,]
    linpr <- lb0 + Z[,1] + (lb1 + Z[,2])*Xdat$x
    p     <- linkinv(linpr)
    nsize <- 25
    Dead  <- rbinom(nrow(Xdat),nsize,p)
    Alive <- nsize - Dead
    x0    <- (linkfun(0.99) - beta0)/beta1
    Xdat$Dead  <- Dead
    Xdat$Alive <- Alive
    attr(Xdat,"trueLD99") <- x0
    return(Xdat)
}

set.seed(42)
X <- artSim()

## 1. Fit with default nlminb control parameters stops too early
system.time(fit1 <- glmmTMB(Dead/(Alive+Dead) ~ (Trt + 0)/x + (x | Rep),
                            weights = Alive+Dead,
                            family=binomial(link="cloglog"),data=X))
logLik(fit1)
range(fit1$sdr$gradient) ## Not converged

## 2. We can fix by increasing the number of interations
con <- glmmTMBControl( optCtrl = list(iter.max=1e3, eval.max=1e3) )
system.time(fit2 <- glmmTMB(Dead/(Alive+Dead) ~ (Trt + 0)/x + (x | Rep),
                            weights = Alive+Dead,
                            family=binomial(link="cloglog"),data=X,
                            control=con))
logLik(fit2)
range(fit2$sdr$gradient) ## Probably OK

## 3. OR we can set profile=TRUE
con <- glmmTMBControl( profile = TRUE )
system.time(fit3 <- glmmTMB(Dead/(Alive+Dead) ~ (Trt + 0)/x + (x | Rep),
                            weights = Alive+Dead,
                            family=binomial(link="cloglog"),data=X,
                            control=con))
logLik(fit3)
range(fit3$sdr$gradient) ## OK
logLik(fit3) > logLik(fit2)

## Compare estimates of fit2 versus fit3
diff <- (summary(fit3)$coefficients$cond - summary(fit2)$coefficients$cond)[,1:2]
range(diff[,"Estimate"])
range(diff[,"Std. Error"])
