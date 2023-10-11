## https://github.com/glmmTMB/glmmTMB/issues/928
## issues with smooths

devtools::load_all("glmmTMB")

## handle NA values properly (use `fr` rather than `mf$data` to adjust smooth)
library(glmmTMB)
data(sleepstudy, package = "lme4")
ss <- sleepstudy
ss[3, "Days"] <- NA
ss$d <- ss$Days+1
m1 <- glmmTMB(Reaction ~ s(Days), sleepstudy)
m2 <- update(m1, data = ss)


m3 <- update(m1, . ~ s(log(Days+1)))
## Error in names(dat) <- object$term : 
##   'names' attribute [1] must be the same length as the vector [0]

m4 <- update(m1, . ~ s(log(d)), data = ss) ## OK

## OK
m5 <- update(m1, . ~ s(log(d), k = 4), data = ss)

f <- function(x, ...) glmmTMB:::findbars_x(x, default.special = NULL, target = "s", ...)
f(Reaction ~ s(log(Days)))
## OK: s(log(Days))
f(Reaction ~ s(log(Days+1)), debug = TRUE)
f(Reaction ~ s(log(d), k = 4))  ## OK


## this does work though ...
f(y ~ a + log(b) + s(x, bs = "tp") + s(y, bs = "gp"))
f(y ~ a + log(b) + s(x, k = 5))  ## works

