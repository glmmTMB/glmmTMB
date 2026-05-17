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
#' @param lambda2 a single numeric value for parameter \code{lambda2} with \eqn{0 \le lambda2 < 1}.
#'                When \code{lambda2=0}, the generalized Poisson distribution
#'                reduces to the Poisson distribution
#' @param log logical; if \code{TRUE}, the log-density is returned
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

#' @rdname dgenpois
#' @export
pgenpois <- function(q, lambda1, lambda2) {
    n <- max(length(q), length(lambda1), length(lambda2))
    q       <- rep_len(q,       n)
    lambda1 <- rep_len(lambda1, n)
    lambda2 <- rep_len(lambda2, n)

    out <- numeric(n)

    for (j in seq_len(n)) {
        qq <- floor(q[j]); l1 <- lambda1[j]; l2 <- lambda2[j]

        if (!is.finite(qq) || qq < 0) { out[j] <- 0; next }
        if (!is.finite(l1) || l1 <= 0 || !is.finite(l2)) { out[j] <- NA_real_; next }

        probs <- dgenpois(0:qq, l1, l2, log = FALSE)
        out[j] <- min(max(sum(probs, na.rm = TRUE), 0), 1)
    }
    out
}


#' @rdname dgenpois
#' @export
rgenpois <-function(n, lambda1, lambda2)
{
  lambda1 <- rep_len(lambda1, n)
  lambda2 <- rep_len(lambda2, n)

  random_genpois <- numeric(n)
  random_unif <- runif(n)
  for (i in 1:n)
  {
    temp_random_genpois <- 0
    kum <- dgenpois(0, lambda1 = lambda1[i], lambda2 = lambda2[i])
    while(random_unif[i] > kum)
    {
      temp_random_genpois <- temp_random_genpois + 1
      kum <- kum + dgenpois(temp_random_genpois, lambda1 = lambda1[i], lambda2 = lambda2[i])
    }
    random_genpois[i] <- temp_random_genpois
  }
  return(random_genpois)
}

## CDF wrappers with (q, mu, ...) signature for use in dunnsmyth_resids switch() blocks.

## pnbinom using the mu parameterization (base R defaults to size/prob).
pnbinom0 <- function(q, mu, ...) pnbinom(q, mu = mu, ...)

## Generalized Poisson CDF: converts glmmTMB's (mu, phi) to (lambda1, lambda2).
## phi here is sigma()^2 as returned by predict(type="disp") for a genpois model.
pgenpois_mu <- function(q, mu, phi, eps = 0.001) {
    phi_val <- sqrt(phi)
    alpha   <- pmin(pmax(1 - 1/phi_val, eps), 1-eps)
    lambda1 <- pmax(mu * (1 - alpha), 1e-10)
    lambda2 <- alpha
    pgenpois(q, lambda1, lambda2)
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
