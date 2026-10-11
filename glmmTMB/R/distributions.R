## implementations taken from HMMpa package, Foraita Ronja and Vitali Winkowski:
## https://github.com/bips-hb/HMMpa

#' The Generalized Poisson Distribution
#'
#' Density, cumulative density, random-deviate functions for the generalized Poisson distribution.
#'
#' @param x a vector of (non-negative integer) quantiles
#' @param q a numeric vector of quantiles
#' @param n number of random deviates to draw
#' @param lambda1 a single numeric value for parameter \code{lambda1} with \eqn{lambda1 > 0}
#' @param lambda2 a single numeric value for parameter \code{lambda2} with \eqn{-1 \le lambda2 < 1}.
#'                When \code{lambda2=0}, the generalized Poisson distribution
#'                reduces to the Poisson distribution
#' @param log logical; if \code{TRUE}, the log-density is returned
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are
#'   \eqn{P(X \le q)}.
#' @param log.p logical; if \code{TRUE}, probabilities are returned on the
#'   log scale.
#'
#' @details
#' The generalized Poisson distribution has the density
#' \deqn{ p(x) = \lambda_1 (\lambda_1 + \lambda_2 \cdot x)^{x-1}
#'   \frac{ \exp(-\lambda_1-\lambda_2 \cdot x) )}{x!}}{%
#'         p(x) = lambda1 (lambda1 + lambda2 x)^(x-1)  exp(-lambda1-lambda2 x) )/x!}
#'   for \eqn{x = 0,1,2,\ldots},
#'   with \eqn{\mbox{E}(X)=
#'   \frac{\lambda_1}{1-\lambda_2}}{E(x)=lambda1/(1-lambda2)} and variance
#'   \eqn{\mbox{var}(X)=\frac{\lambda_1}{(1-\lambda_2)^3}}{var(x)=lambda1/(1-lambda2)^3}.
#'
#' @references Joe, H., Zhu, R. (2005). Generalized Poisson distribution: the property of
#' mixture of Poisson and comparison with negative binomial distribution.
#' Biometrical Journal \bold{47}(2):219--229.
#'
#' @author Based on Joe and Zhu (2005). Original implementation by  Vitali Witowski (2013).
#'
#' @seealso \link{Distributions} for other standard distributions,
#' including \code{\link{dpois}} for the Poisson distribution.
#'
#' @keywords distribution
#'
#' @return
#' \code{dgenpois} gives the density, \code{pgenpois} gives the cumulative
#' distribution function, and \code{rgenpois} generates random deviates.
#'
#' @export
#'
#' @examples
#' dgenpois(x = seq(0, 20), lambda1 = 10, lambda2 = 0.5)
#' pgenpois(q = 5, lambda1 = 10, lambda2 = 0.5)
#' set.seed(101)
#' hist(rgenpois(n = 1000, lambda1 = 10, lambda2 = 0.5))
dgenpois <- function(x, lambda1, lambda2, log = FALSE)
{
    n       <- max(length(x), length(lambda1), length(lambda2))
    x       <- rep_len(x,       n)
    lambda1 <- rep_len(lambda1, n)
    lambda2 <- rep_len(lambda2, n)

    out <- rep(NA_real_, n)
    valid_x <- is.finite(x) & x >= 0 & x == floor(x)
    bad_x <- !valid_x
    out[bad_x] <- if (log) -Inf else 0

    valid_par <- is.finite(lambda1) & lambda1 > 0 & is.finite(lambda2)
    bad_par <- valid_x & !valid_par
    out[bad_par] <- NA_real_

    term <- lambda1 + x * lambda2
    valid_term <- is.finite(term) & term > 0
    bad_term <- valid_x & valid_par & !valid_term
    out[bad_term] <- if (log) -Inf else 0

    ok <- valid_x & valid_par & valid_term
    if (any(ok)) {
        logpmf <- log(lambda1[ok]) + (x[ok] - 1) * log(term[ok]) - term[ok] - lgamma(x[ok] + 1)
        if (log) {
            out[ok] <- logpmf
        } else {
            vals <- exp(logpmf)
            vals[!is.finite(vals)] <- 0
            out[ok] <- vals
        }
    }
    out
}

#' @rdname dgenpois
#' @export
pgenpois <- function(q, lambda1, lambda2, lower.tail = TRUE, log.p = FALSE) {
    n <- max(length(q), length(lambda1), length(lambda2))
    q       <- rep_len(q,       n)
    lambda1 <- rep_len(lambda1, n)
    lambda2 <- rep_len(lambda2, n)

    out <- rep(NA_real_, n)
    q <- floor(q)
    bad_q <- !is.finite(q) | q < 0
    out[bad_q] <- 0

    valid_par <- is.finite(lambda1) & lambda1 > 0 & is.finite(lambda2)
    valid <- !bad_q & valid_par
    out[!bad_q & !valid_par] <- NA_real_

    if (any(valid)) {
      ## hash, find unique values
        key <- paste(lambda1[valid], lambda2[valid], sep = "\t")
        split_idx <- split(which(valid), key)
        for (idx in split_idx) {
            qq <- q[idx]
            max_q <- max(qq)
            probs <- dgenpois(0:max_q, lambda1[idx[1]], lambda2[idx[1]], log = FALSE)
            cdf <- pmin(pmax(cumsum(probs), 0), 1)
            out[idx] <- cdf[qq + 1]
        }
    }
    if (!lower.tail) out <- 1 - out
    if (log.p)       out <- log(out)
    out
}


#' @rdname dgenpois
#' @export
rgenpois <- function(n, lambda1, lambda2) {
    .Call("rgenpois_R",
          as.integer(n),
          as.double(rep_len(lambda1, n)),
          as.double(rep_len(lambda2, n)),
          PACKAGE = "glmmTMB")
}

## CDF wrappers with (q, mu, ...) signature for use in dunnsmyth_resids switch() blocks.

## pnbinom using the mu parameterization (base R defaults to size/prob).
pnbinom0 <- function(q, mu, ...) pnbinom(q, mu = mu, ...)

#' Generalized Poisson CDF: converts glmmTMB's (mu, phi) to (lambda1, lambda2).
#' phi here is sigma()^2 as returned by predict(type="disp") for a genpois model.
#' lambda2 is less than -1 for phi < 0.25. 
#' @param q vector of quantiles
#' @param mu vector of mean values
#' @param phi vector of dispersion parameters
#' @param oor_na return NA values for 'out of range' lambda values (< 0.25)
#' @param oor_clamp clamp out of range values to the appropriate bounds?
#' @param oor_warn_clamp warn when clamping out-of-range values?
#' @param eps epsilon value for clamping
#' @return vector of CDF values
#' @noRd
pgenpois_mu <- function(q, mu, phi, oor_na = TRUE,
                        oor_clamp = FALSE, oor_warn_clamp = TRUE,
                        eps = 1e-6) {
    phi_val <- sqrt(phi)
    alpha   <- 1 - 1/phi_val
    lambda1 <- mu * (1 - alpha)
    lambda2 <- alpha

    invalid <- !is.finite(mu) | !is.finite(phi) | !is.finite(alpha) |
      !is.finite(lambda1)

    if (any(invalid, na.rm = TRUE)) {
        warning(
          "genpois residual CDF parameters are non-finite; ",
          "returning NA for those entries"
        )
    }

    out_of_range <- lambda1 <= 0 | abs(lambda2) >= 1
    if (oor_clamp && any(out_of_range, na.rm = TRUE)) {
      if (oor_warn_clamp) warning("genpois residual CDF parameters out of range: clamping")
      lambda1 <- pmax(eps, lambda1)
      lambda2 <- pmin(1-eps, pmax(lambda2, -(1-eps)))
    }

    if (oor_na && any(out_of_range, na.rm = TRUE)) {
        warning("some genpois parameters out of range: returning NAs")
    }
    out <- rep(NA_real_, max(length(q), length(lambda1), length(lambda2)))
    if (any(!invalid, na.rm = TRUE)) {
        skip <- invalid | (oor_na & out_of_range)
        if (any(!skip)) {
            out[!skip] <- pgenpois(q[!skip], lambda1[!skip], lambda2[!skip])
        }
    }
    out
}

#' The Bell Distribution
#'
#' Density and cumulative distribution function for the Bell distribution.
#'
#' @param x vector of non-negative integer quantiles.
#' @param q vector of quantiles.
#' @param theta vector of positive Bell parameters. Related to the mean by
#'   \eqn{\mu = \theta e^{\theta}}.
#' @param mu vector of positive means. Related to the Bell parameter by
#'   \eqn{\theta = W(\mu)}, where \eqn{W} is the Lambert W function.
#' @param log logical; if \code{TRUE} the log-density is returned.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are
#'   \eqn{P(X \le q)}.
#' @param log.p logical; if \code{TRUE} probabilities are returned on the
#'   log scale.
#'
#' @details
#' The Bell distribution (Castellares et al. 2018) has probability mass function
#' \deqn{P(X = x) = \frac{e^{1 - e^{\theta}} \theta^{x} B_{x}}{x!}}
#' for \eqn{x = 0, 1, 2, \ldots} and \eqn{\theta > 0}, where \eqn{B_x} is the
#' \eqn{x}-th Bell number. The mean is \eqn{\mu = \theta e^{\theta}}.
#'
#' \code{dbell} is parameterized by \code{theta}; \code{pbell} is
#' parameterized by the mean \code{mu} (as used in \code{\link{bell}} models).
#' Both use internal C-level implementations with no external package
#' dependencies.
#'
#' @references Castellares, F., Ferrari, S. L. P., Lemonte, A. J. (2018).
#' On the Bell distribution and its associated regression model for count data.
#' \emph{Applied Mathematical Modelling} \bold{56}:172--185.
#' \doi{10.1016/j.apm.2017.12.014}
#'
#' @seealso \code{\link{bell}} for the Bell family in \pkg{glmmTMB};
#'   \link{Distributions} for other standard distributions.
#'
#' @keywords distribution
#'
#' @return
#' \code{dbell} gives the density; \code{pbell} gives the cumulative
#' distribution function.
#'
#' @examples
#' dbell(0:5, theta = 1)
#' pbell(0:5, mu = exp(1))  ## theta = 1 implies mu = e
#'
#' @export
dbell <- function(x, theta, log = FALSE) {
    .Call("dbell_R", as.double(x), rep_len(as.double(theta), length(x)),
          as.integer(log), PACKAGE = "glmmTMB")
}

#' The Conway-Maxwell-Binomial Distribution
#'
#' Density and cumulative distribution function for the
#' Conway-Maxwell-Binomial (CMB) distribution, parameterized by the mean.
#'
#' @param x vector of quantiles (non-negative integers \eqn{\le} \code{size}).
#' @param q vector of quantiles.
#' @param size vector of numbers of trials.
#' @param mu vector of means (of the \emph{count}, i.e. \eqn{0 < \mu <}
#'   \code{size}).
#' @param nu vector of dispersion parameters; \eqn{\nu = 1} recovers the
#'   binomial distribution, \eqn{\nu > 1} gives under-dispersion,
#'   \eqn{0 < \nu < 1} over-dispersion, and \eqn{\nu < 0} super-dispersion
#'   (see \code{allow_negative_nu} in \code{\link{combinomial}}).
#' @param log logical; if \code{TRUE} the log-density is returned.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are
#'   \eqn{P(X \le q)}.
#' @param log.p logical; if \code{TRUE} probabilities are returned on the
#'   log scale.
#'
#' @details
#' The CMB distribution (Shmueli et al. 2005; Kadane 2016) has probability
#' mass function
#' \deqn{P(X = x) = \binom{n}{x}^{\nu} p^{x} (1-p)^{n-x} / Z(p, \nu, n)}
#' for \eqn{x = 0, 1, \ldots, n}, where \eqn{Z} normalizes over the support.
#' Following the mean parameterization used by the \code{\link{combinomial}}
#' family (cf. Huang 2017), \code{mu} is the mean count; the natural
#' parameter \eqn{p} is recovered internally by solving
#' \eqn{E[X \mid n, p, \nu] = \mu}.
#'
#' @references
#' Shmueli, G., Minka, T. P., Kadane, J. B., Borle, S., Boatwright, P. (2005).
#' A useful distribution for fitting discrete data: revival of the
#' Conway-Maxwell-Poisson distribution.
#' \emph{Journal of the Royal Statistical Society C} \bold{54}(1):127--142.
#'
#' Kadane, J. B. (2016). Sums of possibly associated Bernoulli variables:
#' the Conway-Maxwell-Binomial distribution.
#' \emph{Bayesian Analysis} \bold{11}(2):403--420.
#'
#' @seealso \code{\link{combinomial}} for the CMB family in \pkg{glmmTMB};
#'   \link{Distributions} for other standard distributions, including
#'   \code{\link{dbinom}} for the binomial distribution.
#'
#' @keywords distribution
#'
#' @return
#' \code{dcombinom} gives the density; \code{pcombinom} gives the cumulative
#' distribution function.
#'
#' @examples
#' ## nu = 1 recovers the binomial
#' all.equal(dcombinom(0:6, size = 6, mu = 3, nu = 1),
#'           dbinom(0:6, size = 6, prob = 0.5))
#' pcombinom(0:6, size = 6, mu = 3, nu = 0.5)
#'
#' @export
dcombinom <- function(x, size, mu, nu, log = FALSE) {
    n <- max(length(x), length(size), length(mu), length(nu))
    .Call("dcombinom_R",
          as.double(rep_len(x, n)),
          as.double(rep_len(mu, n)),
          as.double(rep_len(nu, n)),
          as.integer(rep_len(size, n)),
          as.logical(log),
          PACKAGE = "glmmTMB")
}

#' @rdname dcombinom
#' @export
pcombinom <- function(q, size, mu, nu, lower.tail = TRUE, log.p = FALSE) {
    n <- max(length(q), length(size), length(mu), length(nu))
    out <- .Call("pcombinom_R",
                 as.double(rep_len(q, n)),
                 as.double(rep_len(mu, n)),
                 as.double(rep_len(nu, n)),
                 as.integer(rep_len(size, n)),
                 PACKAGE = "glmmTMB")
    if (!lower.tail) out <- 1 - out
    if (log.p)       out <- log(out)
    out
}

#' CMB CDF wrapper for dunnsmyth_resids: converts the proportion-scale
#' (mu, phi) that residuals.glmmTMB works with to the count-scale
#' (mu, nu, size) parameterization of pcombinom. phi is nu as returned
#' by predict(type = "disp") (under either dispersion link).
#' @param q vector of quantiles (count scale)
#' @param mu vector of mean proportions
#' @param phi vector of dispersion (nu) values
#' @param size vector of numbers of trials
#' @return vector of CDF values
#' @noRd
pcombinom_mu <- function(q, mu, phi, size) {
    pcombinom(q, size = size, mu = mu * size, nu = phi)
}

#' @rdname dbell
#' @export
## theta = LambertW(mu); P(X <= q) accumulated from the PMF over 0:floor(q).
pbell <- function(q, mu, lower.tail = TRUE, log.p = FALSE) {
    n     <- max(length(q), length(mu))
    q     <- floor(rep_len(q,  n))
    mu    <- rep_len(mu, n)
    theta <- .Call("lambertW_R", as.double(mu), PACKAGE = "glmmTMB")
    out   <- numeric(n)
    for (j in seq_len(n)) {
        qq <- q[j]
        if (!is.finite(qq) || qq < 0) { out[j] <- 0; next }
        xs     <- as.double(0:qq)
        probs  <- .Call("dbell_R", xs, rep(theta[j], length(xs)), 0L,
                        PACKAGE = "glmmTMB")
        out[j] <- min(max(sum(probs, na.rm = TRUE), 0), 1)
    }
    if (!lower.tail) out <- 1 - out
    if (log.p)       out <- log(out)
    out
}
