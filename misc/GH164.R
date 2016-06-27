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
glmmTMB(y~1,data=d,family=list(family="nbinom2",link="log"))
glmmTMB(y2~x+(1|f),data=d,family=list(family="nbinom2",link="log"))
