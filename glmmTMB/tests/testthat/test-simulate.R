stopifnot(require("testthat"),
          require("glmmTMB"))

set.seed(101)
test_that("simulate", {
    sm2 <<- rowMeans(do.call(cbind, simulate(fm2, 10)))
    sm2P <<- rowMeans(do.call(cbind, simulate(fm2P, 10)))
    sm2G <<- rowMeans(do.call(cbind, simulate(fm2G, 10)))
    sm2NB <<- rowMeans(do.call(cbind, simulate(fm2NB, 10)))
    expect_equal(sm2, sleepstudy$Reaction, tol=20)
	expect_equal(sm2P, sleepstudy$Reaction, tol=20)
	expect_equal(sm2G, sleepstudy$Reaction, tol=20)
	expect_equal(sm2NB, sleepstudy$Reaction, tol=20)
})

test_that("binomial", {
    s1 <- simulate(f1b, 5, seed=1)
    s2 <- simulate(f2b, 5, seed=1)
    s3 <- simulate(f3b, 5, seed=1)
    expect_equal(max(abs(as.matrix(s1) - as.matrix(s2))), 0)
    expect_equal(max(abs(as.matrix(s1) - as.matrix(s3))), 0)
})

test_that("simulate t distrib", {
    data("sleepstudy", package = "lme4")
    sleepstudy <- within(sleepstudy,
                         y <- drop(scale(Reaction, center = min(Reaction),
                                         scale=diff(range(Reaction)))))
    m <- glmmTMB(y ~ Days + (1 | Subject), data = sleepstudy, family = t_family)
    ## GH 1024; sim actually implemented, get different values for each run ...
    set.seed(101)
    expect_equal(
        head(simulate(m, nsim = 1)[[1]], 1),
        0.105283057216147)
    expect_equal(
        head(simulate(m, nsim = 1)[[1]], 1),
        0.342968855368643)
    expect_equal(
        head(simulate(m, nsim = 1, seed = 101)[[1]], 1),
        0.105283057216147)
})
                 
