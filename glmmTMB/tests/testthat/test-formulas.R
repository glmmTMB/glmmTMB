stopifnot(require("testthat"),
          require("glmmTMB"))

context("formula parsing")

nrt <- function(x) length(x$reTrmFormulas)

test_that("basic splitForm", {
  expect_equal(nrt(splitForm(y~(x+q))),0) ## reTrms part should be empty
  sf1 <- splitForm(y~(x+q)+(1|f))
  sf2 <- splitForm(y~(x+q)+us(1|f))
  sf3 <- splitForm(y~(x+q)+diag(1|f))
  sf4 <- splitForm(~x+y+(f|g)+cs(1|g))
  expect_equal(nrt(sf1),1)
  expect_equal(sf1$reTrmFormulas,list(quote(1|f)))
  expect_equal(sf1,sf2) 
  expect_equal(sf3$reTrmClasses,"diag")
  expect_equal(sf4$reTrmClasses,c("us","cs"))
})

test_that("slash terms", {
  sf5 <- splitForm(~x+y+(1|f/g)) 
  sf6 <- splitForm(~x+y+(1|f/g/h))
  sf7 <- splitForm(~x+y+(1|(f/g)/h)) 
  expect_equal(sf5$reTrmClasses,rep("us",2))
  expect_equal(sf6$reTrmClasses,rep("us",3))
  expect_equal(sf6,sf7)
})

test_that("grpvar terms", {
  sf8 <- splitForm(~x+y+(1|f*g)) 
  sf9 <- splitForm(~x+y+(1|f+g+h))
  expect_equal(sf8$reTrmClasses,rep("us",3))
  expect_equal(sf8$reTrmFormula,list(quote(1|f),quote(1|g),quote(1|f:g)))
  expect_equal(sf9$reTrmClasses,rep("us",3))
  expect_equal(sf9$reTrmFormula,list(quote(1|f),quote(1|g),quote(1|h)))
})


test_that("noSpecial", {
    ## handle parentheses in formulas: GH #174
    ff <- y~1+(((us(1|f))))
    expect_equal(noSpecials(ff,delete=FALSE),y~1+(1|f))
    expect_equal(noSpecials(ff),y~1)
    ## 'naked' special - left alone: GH #261
    ff2 <- y ~ us
    expect_equal(noSpecials(ff2),ff2)

})
