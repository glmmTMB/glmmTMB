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

  if (length(x) < max(length(lambda1), length(lambda2)))
  {
  	x <- c(rep(x, times = max(length(lambda1), length(lambda2))))
  }
  if (length(lambda1) < max(length(x), length(lambda2)))
  {
  	lambda1 <- c(rep(lambda1, times = max(length(x), length(lambda2))))
  }
  if (length(lambda2) < max(length(x), length(lambda1)))
  {
  	lambda2 <- c(rep(lambda2, times = max(length(x), length(lambda1))))
  }
  a <- NULL
  for (j in 1:max(length(x), length(lambda1), length(lambda2)))
  {
  	if(x[j] < 2)
  {
  	b = (lambda1[j] * (lambda1[j] + x[j] * lambda2[j])^(x[j] - 1) * (exp(-(
    lambda1[j] + x[j] * lambda2[j])))) / (factorial(x[j]))
  } else  {
  	  f1 <- (lambda1[j] + x[j] * lambda2[j])
      g1 <- exp(-lambda2[j])
      e <- 1
      for (i in 2:x[j])
      {
      	e <- e * ((f1 * g1) / i)
      }
      d1 <- lambda1[j] * exp(-lambda1[j]) * g1
      b <- e * d1
      if (is.na(b))
      {
      	b <- 4.940656e-324
      }
      if(b == Inf)
      {
      	b <- 4.940656e-324
      }
      if(b == -Inf)
      {
      	b <- 4.940656e-324
      }
    }
    a <- c(a,b)
  }
return(a)
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

## Bell CDF: converts mu to the Bell theta parameter via the Lambert W function.
## Requires the 'bellreg' and 'LambertW' packages (listed in Suggests).
pbell_mu <- function(q, mu) {
    if (!requireNamespace("bellreg",  quietly = TRUE) ||
        !requireNamespace("LambertW", quietly = TRUE)) {
        stop("packages 'bellreg' and 'LambertW' are required for Bell family residuals")
    }
    theta <- LambertW::W(mu)
    mapply(function(qq, th) {
        if (qq < 0) return(0)
        val <- suppressWarnings(bellreg::pbell(qq, th))
        if (!is.finite(val)) NA_real_ else val
    }, q, theta)
}
