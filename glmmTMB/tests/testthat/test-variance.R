stopifnot(require("testthat"),
          require("glmmTMB"))

## Check family variance functions (via full_variance()) against the
## variance of simulated responses, and Pearson residuals against their
## expected mean square.
## Each case gives the internal parameters for simulate_new() (eta, betadisp,
## psi) and the same quantities on the scales that full_variance() takes
## (phi as from predict(type = "disp"), shape as from family_params()).
## size: number of trials for binomial-type families.
## zi: zero-inflation probability (cases with zi are zero-inflated, or
##  hurdle models for truncated families).
var_cases <- list(
    list(family = gaussian(), eta = 1, betadisp = log(2), phi = 2),
    list(family = Gamma(link = "log"), eta = log(3), betadisp = log(2),
         phi = 1/sqrt(2)),
    list(family = ziGamma(link = "log"), eta = log(3), betadisp = log(2),
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

## zero-inflated and hurdle models
var_cases <- c(var_cases, list(
    list(family = gaussian(), eta = 1, betadisp = log(2), phi = 2, zi = 0.3),
    list(family = ziGamma(link = "log"), eta = log(3), betadisp = log(2),
         phi = 1/sqrt(2), zi = 0.3),
    list(family = poisson(), eta = log(3), phi = 1, zi = 0.3),
    list(family = nbinom2(), eta = log(3), betadisp = log(2), phi = 2,
         zi = 0.3),
    list(family = truncated_poisson(), eta = log(3), zi = 0.3),
    list(family = truncated_nbinom2(), eta = log(3), betadisp = log(2),
         phi = 2, zi = 0.3)
))

if (requireNamespace("gsl", quietly = TRUE)) {
    var_cases <- c(var_cases, list(list(family = bell(), eta = log(2))))
}

case_name <- function(x) {
    nm <- x$family$family
    if (nm == "compois") nm <- sprintf("%s (phi = %g)", nm, x$phi)
    if (!is.null(x$zi)) nm <- paste(nm, "(zero-inflated)")
    nm
}

sim_case <- function(x, newdata, beta) {
    np <- list(beta = beta)
    if (!is.null(x$betadisp)) np$betadisp <- x$betadisp
    if (!is.null(x$psi)) np$psi <- x$psi
    if (!is.null(x$zi)) np$betazi <- qlogis(x$zi)
    args <- list(~ 1 + z, family = x$family, newdata = newdata,
                 newparams = np, seed = 101)
    if (!is.null(x$zi)) args$ziformula <- ~ 1
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
        ## the hurdle mean square needs the truncated mean; Pearson only
        if (!is.null(x$zi) && grepl("^truncated", x$family$family)) next
        y <- sim_case(x, data.frame(id = seq_len(1e5)), beta = x$eta)
        mu <- x$family$linkinv(x$eta)
        v <- glmmTMB:::full_variance(x$family, mu = mu,
                                     phi = x$phi, shape = x$shape,
                                     size = x$size)
        if (!is.null(x$zi)) v <- v * (1 - x$zi) + mu^2 * x$zi * (1 - x$zi)
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
        if (!is.null(x$zi)) fargs$ziformula <- ~ 1
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
                        ##  (with or without zero-inflation)
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

## independent references: the natural parameter is found with uniroot()
## and the moments are summed over the support
test_that("compois variance matches direct summation", {
    ## P(y) propto lambda^y / (y!)^nu; glmmTMB's phi is 1/nu
    cmp_var <- function(mu, phi, ymax = 1000) {
        nu <- 1/phi
        y <- 0:ymax
        moments <- function(loglambda) {
            lw <- y * loglambda - nu * lfactorial(y)
            w <- exp(lw - max(lw))
            w <- w/sum(w)
            m <- sum(y * w)
            c(mean = m, var = sum((y - m)^2 * w))
        }
        loglambda <- uniroot(function(ll) moments(ll)[["mean"]] - mu,
                             c(-20, 50), tol = 1e-12)$root
        moments(loglambda)[["var"]]
    }
    pars <- expand.grid(mu = c(0.5, 3, 20), phi = c(0.3, 1, 3))
    ref <- mapply(cmp_var, pars$mu, pars$phi)
    expect_equal(compois()$variance(pars$mu, pars$phi), ref,
                 tolerance = 1e-6)
})

test_that("combinomial variance matches direct summation", {
    ## P(y) propto choose(n, y)^nu p^y (1-p)^(n-y);
    ## variance returned per trial, Var(Y)/size
    cmb_var <- function(prob, nu, size) {
        y <- 0:size
        moments <- function(logitp) {
            lw <- nu * lchoose(size, y) + y * logitp
            w <- exp(lw - max(lw))
            w <- w/sum(w)
            m <- sum(y * w)
            c(mean = m, var = sum((y - m)^2 * w))
        }
        logitp <- uniroot(function(lp) moments(lp)[["mean"]] - prob * size,
                          c(-40, 40), tol = 1e-12)$root
        moments(logitp)[["var"]]/size
    }
    pars <- expand.grid(prob = c(0.1, 0.5, 0.8), nu = c(0.3, 1, 3),
                        size = c(5, 20))
    ref <- mapply(cmb_var, pars$prob, pars$nu, pars$size)
    expect_equal(combinomial()$variance(pars$prob, pars$nu, pars$size), ref,
                 tolerance = 1e-6)
    ## nu = 1 is the binomial
    expect_equal(combinomial()$variance(c(0.1, 0.5), 1, 10),
                 binomial()$variance(c(0.1, 0.5)), tolerance = 1e-8)
})
