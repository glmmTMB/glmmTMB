set.seed(101)
dd <- data.frame(y = rnbinom(100, mu = 5, size = 2))
dd <- dd[dd$y >= 3, , drop=FALSE]

library(glmmTMB)
glmmTMB(y ~ 1, family = list(family = "truncated_nbinom2_2", link = "log"),
        data = dd)
