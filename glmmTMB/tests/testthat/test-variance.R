stopifnot(require("testthat"),
          require("glmmTMB"))

## Check family variance functions (via full_variance()) against the
## variance of simulated responses, and Pearson residuals against their
## expected mean square.
## Each case gives the internal parameters for simulate_new() (eta, betadisp,
## psi) and the same quantities on the scales that full_variance() takes
## (phi as from predict(type = "disp"), shape as from family_params()).
## size: number of trials for binomial-type families.
var_cases <- list(
    list(family = gaussian(), eta = 1, betadisp = log(2), phi = 2),
    list(family = Gamma(link = "log"), eta = log(3), betadisp = log(2),
         phi = 1/sqrt(2)),
    list(family = poisson(), eta = log(3), phi = 1),
    list(family = binomial(), eta = qlogis(0.3), phi = 1, size = 10),
    list(family = nbinom2(), eta = log(3), betadisp = log(2), phi = 2),
    list(family = nbinom1(), eta = log(3), betadisp = log(2), phi = 2),
    list(family = nbinom12(), eta = log(3), betadisp = log(2), phi = 2,
         psi = log(1.5), shape = 1.5),
    list(family = compois(), eta = log(3), betadisp = log(2), phi = 2),
    list(family = compois(), eta = log(3), betadisp = log(0.5), phi = 0.5),
    list(family = genpois(), eta = log(3), betadisp = log(2), phi = 2),
    list(family = truncated_poisson(), eta = log(3)),
    list(family = truncated_nbinom2(), eta = log(3), betadisp = log(2),
         phi = 2),
    list(family = truncated_nbinom1(), eta = log(3), betadisp = log(2),
         phi = 2),
    list(family = truncated_genpois(), eta = log(3), betadisp = log(2),
         phi = 2),
    list(family = truncated_compois(), eta = log(3), betadisp = log(2),
         phi = 2),
    list(family = beta_family(), eta = qlogis(0.3), betadisp = log(5),
         phi = 5),
    list(family = betabinomial(), eta = qlogis(0.3), betadisp = log(5),
         phi = 5, size = 10),
    list(family = combinomial(), eta = qlogis(0.3), betadisp = log(0.5),
         phi = 0.5, size = 10),
    list(family = tweedie(), eta = log(3), betadisp = log(2), phi = 2,
         psi = qlogis(0.5), shape = 1.5),
    list(family = skewnormal(), eta = 10, betadisp = log(2), phi = 2,
         psi = 3, shape = 3),
    list(family = lognormal(), eta = log(3), betadisp = log(2), phi = 2),
    list(family = t_family(), eta = 0, betadisp = log(2), phi = 2,
         psi = log(8), shape = 8)
)

if (requireNamespace("gsl", quietly = TRUE)) {
    var_cases <- c(var_cases, list(list(family = bell(), eta = log(2))))
}

case_name <- function(x) {
    nm <- x$family$family
    if (nm == "compois") nm <- sprintf("%s (phi = %g)", nm, x$phi)
    nm
}

sim_case <- function(x, newdata, beta) {
    np <- list(beta = beta)
    if (!is.null(x$betadisp)) np$betadisp <- x$betadisp
    if (!is.null(x$psi)) np$psi <- x$psi
    args <- list(~ 1 + z, family = x$family, newdata = newdata,
                 newparams = np, seed = 101)
    if (length(beta) == 1) args[[1]] <- ~ 1
    if (!is.null(x$size)) args$weights <- rep(x$size, nrow(newdata))
    y <- do.call(simulate_new, args)[[1]]
    ## binomial-type responses on the proportion scale
    if (!is.null(x$size)) y <- y/x$size
    y
}

test_that("family variance functions match simulated variances", {
    skip_on_cran()
    for (x in var_cases) {
        y <- sim_case(x, data.frame(id = seq_len(1e5)), beta = x$eta)
        v <- glmmTMB:::full_variance(x$family, mu = x$family$linkinv(x$eta),
                                     phi = x$phi, shape = x$shape,
                                     size = x$size)
        sim_v <- var(y)
        ## full_variance() gives the per-trial variance for binomial types
        if (!is.null(x$size)) sim_v <- sim_v * x$size
        expect_equal(sim_v, v, tolerance = 0.03, label = case_name(x))
    }
})

test_that("Pearson residuals have the expected mean square", {
    skip_on_cran()
    n <- 5000
    dd <- data.frame(z = seq(-1, 1, length.out = n))
    for (x in var_cases) {
        dd$y <- sim_case(x, dd, beta = c(x$eta, 0.5))
        fargs <- list(y ~ z, family = x$family, data = dd)
        if (!is.null(x$size)) fargs$weights <- rep(x$size, n)
        m <- suppressWarnings(do.call(glmmTMB, fargs))
        r2 <- mean(suppressMessages(residuals(m, type = "pearson"))^2)
        fam <- x$family$family
        expected <- if (fam == "betabinomial") {
                        ## not scaled by the dispersion factor
                        (sigma(m) + x$size)/(sigma(m) + 1)
                    } else if (identical(names(formals(x$family$variance)), "mu") &&
                               fam != "bell") {
                        ## base-R convention: not scaled by sigma
                        sigma(m)^2
                    } else 1
        expect_equal(r2, expected, tolerance = 0.05, label = case_name(x))
    }
})

test_that("variance functions with extra parameters", {
    ## skew-normal: phi is the standard deviation, for any shape
    expect_equal(skewnormal()$variance(1:3, phi = 2), rep(4, 3))
    ## Student t: infinite variance for df <= 2, undefined for df <= 1
    expect_equal(t_family()$variance(0, phi = 2, shape = c(8, 2, 1)),
                 c(4 * 8/6, Inf, NaN))
    expect_error(glmmTMB:::full_variance(nbinom12(), mu = 1, phi = 1),
                 "needs .*psi")
    ## full beta-binomial variance per trial
    expect_equal(glmmTMB:::full_variance(betabinomial(), mu = 0.3, phi = 5,
                                         size = 10),
                 0.3 * 0.7 * 15/6)
})
