## extend examples?
## deal with b-mapping stuff
## conditional REs?
## expose X-construction helper function?

## comparing

library(glmmTMB)
library(mgcv)
library(scam)

data("Nile")
ndat <- data.frame(time = c(time(Nile)), val = c(Nile))

sm1 <- glmmTMB(val ~ s(time, bs = "tp"), data = ndat,
                REML = TRUE, start = list(theta = 5))
sm1_gcv <- mgcv::gam(val ~ s(time, bs = "tp"), data = ndat,
               method = "REML")
plot(val ~ time, data = ndat)
lines(ndat$time, predict(sm1))
lines(ndat$time, predict(sm1_gcv), col = 2)
if (requireNamespace("scam")) {
  sm1M <- glmmTMB(val ~ s(time, bs = "mpd"), data = ndat,
               REML = TRUE, start = list(theta = 5))
  sm1M_scam <- scam(val ~ s(time, bs = "mpd"), data = ndat)
  
  lines(ndat$time, predict(sm1M), lty = 2, col = 4)
  lines(ndat$time, predict(sm1M_scam), lty = 2, col = 4)
}


## save smooth info -- where?
devtools::load_all("glmmTMB")
data("Nile")
ndat <- data.frame(time = c(time(Nile)), val = c(Nile))

sm1 <- glmmTMB(val ~ s(time, bs = "tp"), data = ndat,
                REML = TRUE, start = list(theta = 5))

lapply(sm1$modelInfo$reTrms,
       \(y) lapply(y$smooth_info, \(x) x$re$beta_ind))

lapply(sm1$modelInfo$reTrms,
       \(y) lapply(y$smooth_info, \(x) x$re$b_ind))


