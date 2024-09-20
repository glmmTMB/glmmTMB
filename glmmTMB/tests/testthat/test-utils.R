stopifnot(require("testthat"),
          require("glmmTMB"))

test_that("addForm", {
    expect_equal(addForm(y~x,~1,~z),y~x+1+z)
    expect_warning(addForm(y~x,z~1),
                   "discarding LHS")
})

test_that("noSpecials", {
    expect_equal(noSpecials(y~1+us(1|f)), y~1)
    expect_equal(noSpecials(y~1+us(1|f),delete=FALSE), y~1+(1|f))
    expect_equal(noSpecials(y~us(1|f)), y ~ 1)
    expect_equal(noSpecials(y~us(1|f), delete=FALSE), y~ (1|f))
    expect_equal(noSpecials(y~us+1), y ~ us + 1)
    expect_equal(noSpecials(~us+1), ~ us + 1)
    expect_equal(noSpecials(~1+x+us(1|f), delete=FALSE), ~ 1 + x + (1|f))
})

test_that("extractForm", {
    expect_equal(extractForm(~a+offset(b),quote(offset)),
                 list(quote(offset(b))))
    expect_equal(extractForm(~c,quote(offset)), NULL)
    expect_equal(extractForm(~a+offset(b)+offset(c),quote(offset)),
                 list(quote(offset(b)),quote(offset(c))))
    expect_equal(extractForm(~offset(x),quote(offset)),
                 list(quote(offset(x))))

})

test_that("get_cor", {
  set.seed(145820L)
  n <- 5L
  S <- cov2cor(crossprod(matrix(rnorm(n * n), n, n)))
  R <- chol(S)
  R[] <- R * rep(1 / diag(R), each = n)
  theta <- R[upper.tri(R)]

  x <- S[lower.tri(S)]
  y1 <- get_cor(theta)
  expect_equal(x, y1)
})

test_that("put_cor", {
    round_trip <- function(C) {
        expect_equal(get_cor(put_cor(C), return_val = "mat"), C)
    }
    set.seed(101)
    for (n in 2:10) {
        for (i in 1:10) {
            round_trip(get_cor(rnorm(n*(n-1)/2), return_val = "mat"))
        }
    }
})

test_that("up2date for models with mapped components", {
    m1_map <- glmmTMB(count ~ mined + (1|site),
                  family=poisson, data=Salamanders,
                  map=list(theta=factor(NA)),
                  start=list(theta=log(10)))
    ## need to save and read back in order to exercise
    ##  reconstruction/retaping part of `up2date()`
    fn <- tempfile()
    saveRDS(m1_map, fn)
    m1_mapr <- readRDS(fn)
    m1_mapu <- up2date(m1_mapr)
    expect_equal(vcov(m1_map), vcov(m1_mapu))
})

