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
