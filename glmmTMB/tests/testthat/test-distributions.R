## Reference values computed from HMMpa::dgenpois / HMMpa::pgenpois
## and bellreg::dbell / bellreg::pbell + LambertW::W

test_that("dgenpois matches HMMpa", {

    xs <- c(0, 1, 5, 10, 20)
    expect_equal(
        dgenpois(xs, lambda1 = 3, lambda2 = 0.4),
        c(0.0497870683678639, 0.100119809880978, 0.10528042186071,
          0.0304214008501773, 0.00125955853442515),
        tolerance = 1e-7
    )
})

test_that("pgenpois matches HMMpa", {

    qs <- c(0, 2, 5, 10)
    expect_equal(
        pgenpois(qs, lambda1 = 3, lambda2 = 0.4),
        c(0.0497870683678639, 0.277420277828986,
          0.637262372751267, 0.916366482952780),
        tolerance = 1e-7
    )
})

test_that("pgenpois lower.tail/log.p options work", {
    q <- c(0, 2, 5, 10)
    p <- pgenpois(q, lambda1 = 3, lambda2 = 0.4)
    expect_equal(pgenpois(q, lambda1 = 3, lambda2 = 0.4, lower.tail = FALSE),
                 1 - p)
    expect_equal(pgenpois(q, lambda1 = 3, lambda2 = 0.4, log.p = TRUE),
                 log(p))
})

test_that("dbell matches bellreg", {
    xs    <- c(0, 1, 2, 5, 10)
    theta <- 1.5
    expect_equal(
        dbell(xs, theta),
        c(0.0307554190699851, 0.0461331286049776, 0.0691996929074664,
          0.10120455087717, 0.056680749963415),
        tolerance = 1e-7
    )
})

test_that("pbell matches bellreg", {
    qs <- c(0, 2, 5, 10)
    mu <- 3.0
    expect_equal(
        pbell(qs, mu),
        c(0.156079344333164, 0.491996053501517,
          0.847338719261640, 0.990350951653520),
        tolerance = 1e-6
    )
})

test_that("dcombinom/pcombinom reduce to the binomial at nu = 1", {
    expect_equal(dcombinom(0:12, size = 12, mu = 6, nu = 1),
                 dbinom(0:12, size = 12, prob = 0.5),
                 tolerance = 1e-12)
    expect_equal(dcombinom(0:10, size = 10, mu = 3, nu = 1),
                 dbinom(0:10, size = 10, prob = 0.3),
                 tolerance = 1e-12)
    expect_equal(pcombinom(0:12, size = 12, mu = 6, nu = 1),
                 pbinom(0:12, size = 12, prob = 0.5),
                 tolerance = 1e-12)
})

test_that("dcombinom/pcombinom are consistent across dispersion regimes", {
    ## over-dispersed, under-dispersed, and super-dispersed (nu < 0)
    for (nu in c(0.3, 0.5, 2, 5, -0.5)) {
        d <- dcombinom(0:15, size = 15, mu = 6, nu = nu)
        expect_equal(sum(d), 1, tolerance = 1e-12)
        expect_true(all(d >= 0))
        ## mean parameterization round-trip
        expect_equal(sum((0:15) * d), 6, tolerance = 1e-8)
        p <- pcombinom(0:15, size = 15, mu = 6, nu = nu)
        expect_equal(p, cumsum(d), tolerance = 1e-12)
    }
})

test_that("dcombinom/pcombinom edge cases", {
    expect_equal(pcombinom(-1, size = 12, mu = 6, nu = 1), 0)
    expect_equal(pcombinom(12, size = 12, mu = 6, nu = 1), 1)
    expect_equal(pcombinom(99, size = 12, mu = 6, nu = 1), 1)
    expect_equal(dcombinom(3.5, size = 12, mu = 6, nu = 1), 0)
    expect_equal(dcombinom(-1, size = 12, mu = 6, nu = 1), 0)
    ## mean outside (0, size) is invalid
    expect_true(is.na(dcombinom(3, size = 12, mu = 0, nu = 1)))
    expect_true(is.na(pcombinom(3, size = 12, mu = 13, nu = 1)))
    ## log/lower.tail options
    p <- pcombinom(c(0, 3, 7), size = 12, mu = 6, nu = 0.5)
    expect_equal(pcombinom(c(0, 3, 7), size = 12, mu = 6, nu = 0.5,
                           lower.tail = FALSE), 1 - p)
    expect_equal(pcombinom(c(0, 3, 7), size = 12, mu = 6, nu = 0.5,
                           log.p = TRUE), log(p))
    expect_equal(dcombinom(0:5, size = 12, mu = 6, nu = 0.5, log = TRUE),
                 log(dcombinom(0:5, size = 12, mu = 6, nu = 0.5)))
})
