library(testthat)
library(glmmTMB)

## test cases: ZI and non-ZI, Poisson and truncated Poisson

set.seed(42)
N <- 100
df <- data.frame(p1 = rpois(N, 1),
                 nb1 = rnbinom(N, mu = 1, size = 1),
                 x = rnorm(N)) |>
    transform(
        ## zero-inflated versions
        zp1 = p1 * rbinom(N, size = 1, prob = 0.5),
        znb1 = nb1 * rbinom(N, size = 1, prob = 0.5),
        ## truncated versions (NAs will be dropped)
        tp1 = ifelse(p1 == 0, NA, p1),
        tnb1 = ifelse(nb1 == 0, NA, nb1))

testfun <- function(model, response, distrib) {
    zp1 <- predict(model, type="zprob")
    cm1 <- predict(model, type="conditional")
    mu1 <- predict(model, type="response")

    nzvals <- df[[response]]>0

    ## easy: zero-inflation probability == fraction zero
    ## (true for all truncated distributions)
    expect_equal(mean(nzvals),
                 mean(1-zp1),
                 tolerance = 1e-5)

    ## mean of non-zero values == conditional mean
    ## (for all truncated distributions)
    expect_equal(mean(df[[response]][nzvals]),
    mean(cm1),
    tolerance = 1e-2)

    ## compute zero-trunc by hand
    eta <- predict(model, type = "link")
    cm2 <- exp(eta)/(1-distrib(0, exp(eta), sigma(model)))
    expect_equal(cm1, cm2)

    ## mean of all values == response
    ## expect_equal(mean(df[[response]]),
    ##              mean(mu1),
    ##              tolerance = 5e-3)
}

my_dpois <- function(x, lambda, ...) dpois(x, lambda)
my_nb2 <- function(x, mu, size) dnbinom(x, mu = mu, size = size)

f_zp1 <- glmmTMB(zp1 ~ x,
                 zi= ~ 1,
                 family=truncated_poisson(link="log"),
                 data=df)


testfun(f_zp1, "zp1", my_dpois)

f_znb2 <- update(f_zp1, znb1 ~ .,
                 family = truncated_nbinom2)

## hmm, fails for nbinom2??
testfun(f_znb2, "znb1", my_nb2)
## calculate by hand
cm1 <- predict(f_znb2, type="conditional")
eta <- predict(f_znb2, type = "link")
## mean value
mval <- function(mu, size) {
    mu/(1-dnbinom(0, mu = mu, size = size))
}
## this seems to work, so something about the identity is off
all.equal(mval(exp(eta), sigma(f_znb2)), cm1)

f_znb1 <- update(f_zp1, znb1 ~ .,
                 family = truncated_nbinom1)
testfun(f_znb1, "znb1")  ## works.

f_tp1 <- update(f_zp1, tp1 ~ .,
                 ziformula = ~0,
                 family = truncated_poisson)
expect_equal(predict(f_tp1, type = "conditional"),
             predict(f_tp1, type = "response"))
expect_true(all(predict(f_tp1, type = "zprob") == 0))
