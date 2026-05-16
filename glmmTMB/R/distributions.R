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
#' @param lambda2 a single numeric value for parameter \code{lambda2} with \eqn{0 \le lamdba2 < 1}.
#'                When \code{lambda2=0}, the generalized Poisson distribution
#'                reduces to the Poisson distribution
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
#' @references Joe, H., Zhu, R. (2005). Generalized poisson distribution: the property of
#' mixture of poisson and comparison with negative binomial distribution.
#' Biometrical Journal \bold{47}(2):219--229.
#'
#' @author Based on Joe and Zhu (2005). Implementation by  Vitali Witowski (2013).
#'
#' @seealso \code{\link{pgenpois}}, \code{\link{rgenpois}};
#'  \link{Distributions} for other standard distributions,
#'  including \code{\link{dpois}} for the Poisson distribution.
#'
#' @keywords distribution
#'
#' @return
#'  \code{\link{dgenpois}} gives the density of the generalized Poisson distribution.
#' @export
#'
#' @examples
#' dgenpois(x = seq(0,20), lambda1 = 10, lambda2 = 0.5)
#' pgenpois(q = 5, lambda1 = 10, lambda2 = 0.5)
#' hist(rgenpois(n = 1000, lambda1 = 10, lambda2 = 0.5) )
dgenpois <- function(x, lambda1, lambda2)
{
    n       <- max(length(x), length(lambda1), length(lambda2))
    x       <- rep_len(x,       n)
    lambda1 <- rep_len(lambda1, n)
    lambda2 <- rep_len(lambda2, n)
    out <- numeric(n)
    for (j in seq_len(n)) {
        if (x[j] < 2) {
            out[j] <- (lambda1[j] * (lambda1[j] + x[j] * lambda2[j])^(x[j] - 1) *
                       exp(-(lambda1[j] + x[j] * lambda2[j]))) / factorial(x[j])
        } else {
            f1 <- lambda1[j] + x[j] * lambda2[j]
            g1 <- exp(-lambda2[j])
            e  <- 1
            for (i in 2:x[j]) e <- e * (f1 * g1 / i)
            b  <- lambda1[j] * exp(-lambda1[j]) * g1 * e
            out[j] <- if (!is.finite(b) || is.na(b)) .Machine$double.xmin else b
        }
    }
    out
}

#' @rdname dgenpois
pgenpois <- function(q, lambda1, lambda2)
  {
    foo <- 0
    for (i in 0:q)
    {
      foo <- foo + dgenpois(i, lambda1 = lambda1, lambda2 = lambda2)
    }
    return(foo)
  }


#' @rdname dgenpois
rgenpois <-function(n, lambda1, lambda2)
{
  random_genpois <- numeric(n)
  for (i in 1:n)
  {
    temp_random_genpois <- 0
    random_number <- runif(1)
    kum <- dgenpois(0, lambda1 = lambda1, lambda2 = lambda2)
    while(random_number > kum)
    {
      temp_random_genpois <- temp_random_genpois + 1
      kum <- kum + dgenpois(temp_random_genpois, lambda1 = lambda1, lambda2 = lambda2)
    }
    random_genpois[i] <- temp_random_genpois
  }
  return(random_genpois)
}

## Internal generalized Poisson PMF/CDF for Dunn-Smyth residuals.
## Uses log-scale arithmetic for numerical robustness; supports vectorized inputs.
##
## Parameterization: P(Y=y) = lambda1*(lambda1 + lambda2*y)^(y-1)*exp(-(lambda1+lambda2*y))/y!
dgenpois_custom <- function(x, lambda1, lambda2, log = FALSE) {
    n <- max(length(x), length(lambda1), length(lambda2))
    x       <- rep_len(x,       n)
    lambda1 <- rep_len(lambda1, n)
    lambda2 <- rep_len(lambda2, n)

    out <- rep(NA_real_, n)

    for (j in seq_len(n)) {
        xx <- x[j]; l1 <- lambda1[j]; l2 <- lambda2[j]

        if (!is.finite(xx) || xx < 0 || xx != floor(xx)) {
            out[j] <- if (log) -Inf else 0
            next
        }
        if (!is.finite(l1) || l1 <= 0 || !is.finite(l2)) {
            out[j] <- NA_real_
            next
        }
        term <- l1 + xx * l2
        if (!is.finite(term) || term <= 0) {
            out[j] <- if (log) -Inf else 0
            next
        }
        logpmf <- log(l1) + (xx - 1) * log(term) - term - lgamma(xx + 1)
        out[j] <- if (log) logpmf else { val <- exp(logpmf); if (is.finite(val)) val else 0 }
    }
    out
}

pgenpois_custom <- function(q, lambda1, lambda2) {
    n <- max(length(q), length(lambda1), length(lambda2))
    q       <- rep_len(q,       n)
    lambda1 <- rep_len(lambda1, n)
    lambda2 <- rep_len(lambda2, n)

    out <- numeric(n)

    for (j in seq_len(n)) {
        qq <- floor(q[j]); l1 <- lambda1[j]; l2 <- lambda2[j]

        if (!is.finite(qq) || qq < 0) { out[j] <- 0; next }
        if (!is.finite(l1) || l1 <= 0 || !is.finite(l2)) { out[j] <- NA_real_; next }

        probs <- dgenpois_custom(0:qq, l1, l2, log = FALSE)
        out[j] <- min(max(sum(probs, na.rm = TRUE), 0), 1)
    }
    out
}

## CDF wrappers with (q, mu, ...) signature for use in dunnsmyth_resids switch() blocks.

## pnbinom using the mu parameterization (base R defaults to size/prob).
pnbinom0 <- function(q, mu, ...) pnbinom(q, mu = mu, ...)

## Generalized Poisson CDF: converts glmmTMB's (mu, phi) to (lambda1, lambda2).
## phi here is sigma()^2 as returned by predict(type="disp") for a genpois model.
pgenpois_mu <- function(q, mu, phi) {
    phi_val <- sqrt(phi)
    alpha   <- pmin(pmax(1 - 1/phi_val, -0.99), 0.999)
    lambda1 <- pmax(mu * (1 - alpha), 1e-10)
    lambda2 <- alpha
    pgenpois_custom(q, lambda1, lambda2)
}

## Bell CDF using C-level dbell_R and lambertW_R (no external package dependencies).
## theta = LambertW(mu); P(X <= q) is accumulated from the PMF over 0:floor(q).
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
