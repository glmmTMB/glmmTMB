set.seed(101)
dd <- data.frame(y = rnbinom(10000, mu = 5, size = 2))
dd <- dd[dd$y >= 3, , drop=FALSE]

## set up family: still need to work out the variance of the 2-truncated nbinom if we want
## everything to work properly (e.g. Pearson residuals)
## conditional variance of a truncated distribution?
##  could do it by brute force:
##   https://data.princeton.edu/wws509/notes/countmoments
## NOT correct: need to go through this in more detail
##  variance = function(mu, phi, trunc_val = 2) {
##    nb_var <- mu * (1 + mu / phi)
##    trunc_mass <- pnbinom(trunc_val, mu = mu, size = phi)
##    i <- sseq(0, trunc_val)
##    mean_miss <- sum(i * dnbinom(i, mu  = mu, size = phi))
##    var_miss <- sum(i^2 * dnbinom(i, mu  = mu, size = phi))
##    return((nb_var - var_miss + mean_miss^2)/(1-trunc_mass))
## })

ff <- glmmTMB:::family_factory("log", "truncated_nbinom2_2", variance = NULL)

library(glmmTMB)
m1 <- glmmTMB(y ~ 1,
              family = ff,
              ## family = list(family = "truncated_nbinom2_2", link = "log"),
        data = dd)
exp(fixef(m1)$cond) ## 4.95 ~ 5 (true value)
sigma(m1)  ## 2.07 ~ 2 (true value)

