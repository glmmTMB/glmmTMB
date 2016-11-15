set.seed(101)
d <- data.frame(y=rnbinom(1e6,size=1,mu=5))
r <- 3*rnorm(1000)
d <- within(d, {
            f <- factor(rep(1:1000,each=1000))
            x <- runif(1e6)
            eta <- with(d,2+3*x+r[f])
            mu <- exp(eta)
            y2 <- rnbinom(1e6,size=1,mu=mu)
            })

library(glmmTMB)
system.time(m1 <- glmmTMB(y~1,data=d,
                   family=list(family="nbinom2",link="log"),verbose=TRUE))
## 39 seconds
system.time(m2 <- update(m1,y2~x+(1|f)))
## Order of parameters:
## [1] "beta"    "b"       "betazi"  "bzi"     "theta"   "thetazi" "betad"  
## Not matching template order:
## [1] "beta"    "betazi"  "b"       "bzi"     "betad"   "theta"   "thetazi"
## Your parameter list has been re-ordered.
## (Disable this warning with checkParameterOrder=FALSE)
## Constructing atomic invpd
## Constructing atomic invpd
## Constructing atomic D_lgamma
## Constructing atomic invpd
## Constructing atomic D_lgamma
## Error in sparseHessianFun(env, skipFixedEffects = skipFixedEffects) : 
##   Memory allocation fail in function 'MakeADHessObject2'


## Error in sparseHessianFun(env, skipFixedEffects = skipFixedEffects) : 
##   Memory allocation fail in function 'MakeADHessObject2'

## In addition: Warning messages:
## 1: In nlminb(start = par, objective = fn, gradient = gr) :
##   NA/NaN function evaluation
## 2: In he(par) : restarting interrupted promise evaluation
## outer mgc:  NaN 
## Error in nlminb(start = par, objective = fn, gradient = gr) : 
##   gradient function must return a numeric vector of length 4
## Timing stopped at: 35.02 54.44 306.3
## Timing stopped at: 68.77 70.27 392.1

