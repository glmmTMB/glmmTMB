library(lme4)
library(glmmTMB)

## set up example
set.seed(123)
warpbreaks$ID <- sample(1:5, nrow(warpbreaks), replace = TRUE)

system.time(m1 <- glmmTMB::glmmTMB(breaks ~ wool * tension + (1 | ID),
                       family = Gamma(link = "inverse"), data = warpbreaks))

## need to set iter.max/eval.max in glmmTMBControl ...
out1 <- capture.output(
    m1B <- update(m1,
           control = glmmTMBControl(optCtrl = list(trace = 1, iter.max=300, eval.max=400))
           )
)

##> 8.75 sec elapsed
## 5.6 seconds for me

## something very weird going on with 

## convert output to numerical matrix
pars1 <- (out1
    |> stringr::str_replace(":", "")
    |> stringr::str_extract_all("[0-9.e-]+")
    |> lapply(as.numeric)
    |> do.call(what = rbind)
)
colnames(pars1) <- c("iter", "nll",
                     names(fixef(m1)$cond),
                     "logdisp", "logsd")

## plot
pfun <- function(x,y, ...) {
    lines(x, y)
    points(x,y, pch=3, col = "gray")
    points(x[1], y[1], col = "black", pch=16)
    points(tail(x,1), tail(y,1), col = "red", pch=16)
}
par(las=1, bty ="l")

pairs(pars1[,-(1:2)], gap = 0, panel = pfun)

## glmer fit
system.time(
    out2 <- capture.output(
        m2 <- lme4::glmer(breaks ~ wool * tension + (1 | ID), family = Gamma(), data = warpbreaks, verbose = 1,
                          control = glmerControl(optCtrl = list(iprint=1)))
    )
)
## 0.1 seconds

## output conversion etc.
pars2 <- (out2
    |> stringr::str_subset("NM")
    |> stringr::str_extract_all("[0-9.e-]+")
)
pars2 <- pars2[lengths(pars2)==9]
pars2 <- (pars2
    |> lapply(as.numeric)
    |> do.call(what = rbind)
)
colnames(pars2) <- c("iter", "nll",
                     "sd", names(fixef(m2)))

pairs(pars2[,-(1:2)], gap = 0, panel = pfun)

all.equal(c1 <- fixef(m1)$cond, c2 <- fixef(m2))
p1 <- c(c1, exp(getME(m1, "theta")))
p2 <- unlist(getME(m2, c("beta", "theta")))

logLik(m1) - logLik(m2)  ## log-lik difference of 0.5; m1 is better

f1 <- getME(m2, "devfun")
f2 <- m1$obj$fn
try(f1(p1))

## TO DO: match parameters from m2 (need to get sigma and convert to log-dispersion scale),
##  try both in f2 for comparison

## compare trajectorie of fitting

plot(nll-min(nll) + 1e-2 ~ iter, data = as.data.frame(pars1), type ="b",
     log = "y",
     xlim = c(0,331))
with(as.data.frame(pars2), lines(iter, (nll-min(nll))/2 + 1e-2, col = "red", type = "b"))
## 'nll' is actually 2*nll for glmer
## glmer profiles out dispersion parameter, so starts out much better
## 


###
source("https://raw.githubusercontent.com/glmmTMB/glmmTMB/master/misc/allFit.R")
af <- allFit(m1)
## mostly ran out of iterations ...
badfit <- sapply(af, function(x) inherits(x, "try-error"))
af_OK <- af[!badfit]
sapply(af_OK, function(x) x$obj$fn())
sapply(af_OK, \(x) attr(x, "system.time")["elapsed"])
sapply(af_OK, \(x) x$fit$convergence)
