## Demonstrate how to use AR1:
## Try term "ar1(Days | Subject)" on sleepstudy example.

## Get example struct we can use as starting point:
library(glmmTMB)
data(sleepstudy, package="lme4")
if(TRUE){
    ## Original example is 'too easy'. Try with multiple days within
    ## subject and randomly re-order the data. Then run multiple times
    ## and check that 'fit' is un-changed.
    levels(sleepstudy$Subject)[1:2] <- "collect308:309"
    sleepstudy <- sleepstudy[sample(nrow(sleepstudy)),]
}
struct <- glmmTMB(Reaction ~ 1, sleepstudy, debug=TRUE)

## Build term "ar1(Days | Subject)".
## Do the book keeping (Must account for un-ordered and/or duplicate
## time points...). First, sort by group and by time within group):
helperDF <- data.frame(obs.id = 1:nrow(sleepstudy),
                       group  = sleepstudy[[ "Subject" ]],
                       time   = sleepstudy[[ "Days"    ]] )
## Sort by time:
helperDF <- helperDF[order(helperDF$time), ]
## Sort by group:
helperDF <- helperDF[order(helperDF$group), ]
## Factor to lookup random effect:
helperDF <- transform(helperDF, lookup = factor(group:factor(time)) )

## From this, build the terms:
timesByGroup <- split(helperDF$time, helperDF$group)
ar1Term <- function(i){
    times <- unique( timesByGroup[[i]] )
    list(blockReps = 1,
         blockSize = length(times),
         ## Note: 0 ==> use same parms as previous term
         blockNumTheta = if (i==1) 2 else 0,
         blockCode = glmmTMB:::.valid_covstruct["ar1"],
         times=times)
}
terms <- lapply(1:length(timesByGroup), ar1Term)

## Build the Z matrix (and permute back):
Z <- model.matrix( ~ lookup - 1, data=helperDF)
Z <- Z[Matrix::invPerm(helperDF$obs.id), , drop=FALSE]
Z <- as(Z, "dgTMatrix")

## New random effects
b <- rep(0, ncol(Z))

## New theta parameters
theta <- rep(0, 2)

## Modify struct with new term/parameters
struct$data.tmb$Z       <- cbind2(struct$data.tmb$Z, Z)
struct$data.tmb$terms   <- c(struct$data.tmb$terms, terms)
struct$parameters$b     <- c(struct$parameters$b, b)
struct$parameters$theta <- c(struct$parameters$theta, theta)

## Run glmmTMB:
library(TMB)
obj <- with(struct, MakeADFun(data.tmb, parameters, random = "b", 
                                profile = NULL, silent = FALSE, DLL = "glmmTMB"))
optTime <- system.time(fit <- with(obj, nlminb(start = par, 
                                               objective = fn, gradient = gr)))

## Check output
(fit)
obj$report()$corr[[1]] ## Estimate corr ~ 1 Day
obj$report()$sd[[1]]   ## Estimate sd of AR1
plot(obj$env$parList()$b, type="b")

## =============================================
## Replicate this example with formula interface:
sleepstudy$DaysFac <- factor(sleepstudy$Days)
fm <- glmmTMB(Reaction ~ 1 + ar1(DaysFac + 0 | Subject), sleepstudy)
all.equal( fit$par, fm$fit$par )
