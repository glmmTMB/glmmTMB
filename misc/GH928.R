## https://github.com/glmmTMB/glmmTMB/issues/928
## issues with smooths

## devtools::load_all("glmmTMB")

## handle NA values properly (use `fr` rather than `mf$data` to adjust smooth)
library(glmmTMB)
ff <- function(x, ...) glmmTMB:::findbars_x(x, default.special = NULL, target = "s", ...)

data(sleepstudy, package = "lme4")
ss <- sleepstudy
ss[3, "Days"] <- NA
ss$d <- ss$Days+1
m1 <- glmmTMB(Reaction ~ s(Days), sleepstudy)
m2 <- update(m1, data = ss)
m3 <- update(m1, . ~ s(log(Days+1)))
m4 <- update(m1, . ~ s(log(d)), data = ss) ## OK
m5 <- update(m1, . ~ s(log(d), k = 4), data = ss)


d <- data.frame(z = as.vector(volcano),
                x = as.vector(row(volcano)),
                y = as.vector(col(volcano)))
set.seed(1)
d$z <- d$z + rnorm(length(volcano), sd=15)
d <- d[sample(nrow(d), 100), ]

d$pos <- numFactor(d$x, d$y)
d$group <- factor(rep(1, nrow(d)))

## f <- glmmTMB(z ~ 1 + exp(pos + 0 | group), data=d,
##             control = glmmTMBControl(parallel = 10))

## does something but not the right thing?
f2 <- glmmTMB(z ~ 1 + s(x, y), data=d,
             control = glmmTMBControl(parallel = 10))

ff(z ~ 1 + s(x, y, bs = "tp"))
## error
f3 <- glmmTMB(z ~ 1 + s(x, y, bs = "tp"), data=d,
             control = glmmTMBControl(parallel = 10))

form <- z ~ 1 + s(x, y, bs = "tp")
sub_specials(form)

f(Reaction ~ s(log(Days)))
## OK: s(log(Days))
f(Reaction ~ s(log(Days+1)), debug = TRUE)
f(Reaction ~ s(log(d), k = 4))  ## OK


## this does work though ...
f(y ~ a + log(b) + s(x, bs = "tp") + s(y, bs = "gp"))
f(y ~ a + log(b) + s(x, k = 5))  ## works
