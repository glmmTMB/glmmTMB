library(glmmTMB)
devtools::load_all("glmmTMB")
## debug(getReStruc)
## undebug(glmmTMB:::mkTMBStruc)
data("sleepstudy", package = "lme4")

## machinery only works with diag() at the moment ...
g <- glmmTMB(Reaction ~ 1 + diag(1+Days|Subject), sleepstudy)

## for now we will set the simulation behaviour by manually modifying the values inside
## the 'data' object in the TMB object's 'env' (environment) component
## eventually, we want to set up simulate() so it does this for us automatically
## (with a new argument that controls this behaviour term-by-term)
## we also need to make sure simulate_new() can access these controls

## abbreviate for convenience
ee <- g$obj$env
## the default simcode ("random" is set for all terms for a new model; this is back-compatible)
sapply(ee$data$terms, \(x) names(.valid_simcode)[match(x$simCode, .valid_simcode)])

## simulate once
s1 <- simulate(g, seed=101)[[1]]
## save b value (b values may?? get overwritten by zero/newly simulated values if not careful; we can restore them if necessary)
b <- getME(g, "b")

## modify simCode for *all* terms (modifies the object *in place*
set_simcodes <- function(g, val = "zero") {
    ee <- g$obj$env
    for (i in seq_along(ee$data$terms)) ee$data$terms[[i]]$simCode <- .valid_simcode[[val]]
}

set_simcodes(g, "zero")
s2 <- simulate(g, seed=101)[[1]]

## restore original
set_simcodes(g,"random")
s3 <- simulate(g, seed=101)[[1]]
stopifnot(identical(s1, s3))  ## we haven't broken anything

set_simcodes(g,"fix")
s4 <- simulate(g, seed=101)[[1]]

## they're different when we expect them to be, so that's good
head(cbind(random = s1, zero = s2, random2 = s3, fix = s4))

## refit (I may have broken/overwritten something)
g <- glmmTMB(Reaction ~ 1 + diag(1+Days|Subject), sleepstudy)
## back to original
set_simcodes(g,"random")
ss1 <- simulate(g, seed = 101, nsim = 100)

set_simcodes(g,"zero")
ss2 <- simulate(g, seed = 101, nsim = 100)

set_simcodes(g,"fix")
ss3 <- simulate(g, seed = 101, nsim = 100)

par(las=1, bty = "l")
matplot(ss1, type = "p", pch = 1, col = adjustcolor("black", alpha.f = 0.2))
matpoints(ss2, pch = 1, col = adjustcolor("red", alpha.f = 0.2))
matpoints(ss3, pch = 1, col = adjustcolor("blue", alpha.f = 0.2))
legend("bottomleft", col=c("black", "red", "blue"), pch = 1,
       legend = c("random", "zero", "fix"))

## the patterns are as we'd expect ...
