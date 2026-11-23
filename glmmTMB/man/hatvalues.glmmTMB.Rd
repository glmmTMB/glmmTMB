\name{hatvalues.glmmTMB}
\alias{hatvalues.glmmTMB}
\title{Hat values (and other influence measures) for glmmTMB objects}
\usage{
\method{hatvalues}{glmmTMB}(model, diag = TRUE, type = c("full", "conditional"), ...)
}
\arguments{
\item{model}{a \code{glmmTMB} model}

\item{diag}{Get diagonal only?}

\item{type}{\code{"full"}: all top-level parameters (fixed effects,
covariance and dispersion parameters) are treated as functions of the data.
\code{"conditional"}: covariance, dispersion and other non-fixed-effect
top-level parameters are held at their estimates, and the fixed effects
(conditional and zero-inflation) are treated jointly with the random effects,
i.e. as the joint conditional mode given the covariance parameters (the
definition used by \code{lme4}). See Details.}

\item{\dots}{unused: for method compatibility only}
}
\description{
Hat values (and other influence measures) for glmmTMB objects
}
\details{
The hat values are the derivatives of the fitted values (on the response
scale) with respect to the corresponding observations.

The conditional version is much faster and uses much less memory when the
model has many covariance parameters (e.g. reduced-rank (\code{rr}) terms),
since the memory used by the full version scales with (number of
observations) x (number of top-level parameters).

For non-Gaussian responses the fixed-effect estimates are not exactly the
joint conditional mode (glmmTMB maximizes the Laplace-approximated marginal
likelihood), so the conditional hat values are an approximation to the
leverage holding the covariance parameters fixed.

Computing hat values currently requires the \pkg{RTMB} package.
}
\examples{
fm0 <- glmmTMB(weight ~ Time, data=ChickWeight)
fm1 <- lm(weight ~ Time, data = ChickWeight)
all.equal(hatvalues(fm0), unname(hatvalues(fm1)))
fm2 <- glmmTMB(weight ~ Time + diag(Time | Chick), data=ChickWeight)
fm3 <- lme4::lmer(weight ~ Time + (Time || Chick), data = ChickWeight, REML = FALSE)
all.equal(hatvalues(fm2, type = "conditional"), unname(hatvalues(fm3)), tolerance = 1e-5)
}
