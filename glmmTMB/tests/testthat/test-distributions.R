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
