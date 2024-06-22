## see https://github.com/glmmTMB/glmmTMB/issues/928 for discussion

library(glmmTMB)
d <- data.frame(z = as.vector(volcano),
                x = as.vector(row(volcano)),
                y = as.vector(col(volcano)))
set.seed(1)
d$z <- d$z + rnorm(length(volcano), sd=15)
d <- d[sample(nrow(d), 100), ]

d$pos <- numFactor(d$x, d$y)
d$group <- factor(rep(1, nrow(d)))

system.time(
    f <- glmmTMB(z ~ 1 + exp(pos + 0 | group), data=d,
                 control = glmmTMBControl(parallel = 10))
)
predict(f)

f2 <- glmmTMB(z ~ 1 + s(x, y), data=d,
             control = glmmTMBControl(parallel = 10))
predict(f2)

f3 <- glmmTMB(z ~ 1 + s(x, y, bs = "tp"), data=d,
             control = glmmTMBControl(parallel = 10))

identical(predict(f2), predict(f3))

pairs(cbind(z = d$z, exp = predict(f), s = predict(f2)),
      gap = 0)
