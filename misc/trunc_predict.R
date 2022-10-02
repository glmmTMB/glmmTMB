library(testthat)

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

f_zp1 <- glmmTMB(zp1 ~ x,
                 zi= ~ 1,
                 family=truncated_poisson(link="log"),
                 data=df)

zp1 <- predict(f_zp1, type="zprob")
cm1 <- predict(f_zp1, type="conditional")
mu1 <- predict(f_zp1, type="response")

nzvals <- df$zp1>0

## easy: zero-inflation probability == fraction zero
## (true for all truncated distributions)
expect_equal(mean(nzvals),
             mean(1-zp1),
             tolerance = 1e-5)

## mean of non-zero values == conditional mean
## (for all truncated distributions)
expect_equal(mean(df$zp1[nzvals]),
             mean(cm1),
             tolerance = 1e-2)

## mean of all values == response
expect_equal(mean(df$zp1),
             mean(mu1),
             tolerance = 5e-3)

