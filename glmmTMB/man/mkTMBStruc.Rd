% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/glmmTMB.R
\name{mkTMBStruc}
\alias{mkTMBStruc}
\title{Extract info from formulas, reTrms, etc., format for TMB}
\usage{
mkTMBStruc(
  formula,
  ziformula,
  dispformula,
  combForm,
  mf,
  fr,
  yobs,
  respCol,
  weights = NULL,
  contrasts,
  family,
  se = NULL,
  call = NULL,
  verbose = NULL,
  ziPredictCode = "corrected",
  doPredict = 0,
  whichPredict = integer(0),
  aggregate = factor(),
  REML = FALSE,
  start = NULL,
  map = NULL,
  sparseX = NULL,
  control = glmmTMBControl(),
  old_smooths = NULL,
  priors = NULL
)
}
\arguments{
\item{formula}{combined fixed and random effects formula, following lme4 syntax.}

\item{ziformula}{a \emph{one-sided} (i.e., no response variable) formula for zero-inflation combining fixed and random effects: the default \code{~0} specifies no zero-inflation. Specifying \code{~.} sets the zero-inflation formula identical to the right-hand side of \code{formula} (i.e., the conditional effects formula); terms can also be added or subtracted. \strong{When using \code{~.} as the zero-inflation formula in models where the conditional effects formula contains an offset term, the offset term will automatically be dropped}. The zero-inflation model uses a logit link.}

\item{dispformula}{a \emph{one-sided} formula for dispersion combining fixed and random effects: the default \code{~1} specifies the standard dispersion given any family. The argument is ignored for families that do not have a dispersion parameter. For an explanation of the dispersion parameter for each family, see \code{\link{sigma}}. The dispersion model uses a log link. In Gaussian mixed models, \code{dispformula=~0} fixes the residual variance to be 0 (actually a small non-zero value), forcing variance into the random effects. The precise value can be controlled via \code{control=glmmTMBControl(zero_dispval=...)}; the default value is \code{sqrt(.Machine$double.eps)}.}

\item{combForm}{combined formula}

\item{mf}{call to model frame}

\item{fr}{model frame}

\item{yobs}{observed y}

\item{respCol}{response column}

\item{weights}{model weights (for binomial-type models, used as size/number of trials)}

\item{contrasts}{an optional list, e.g., \code{list(fac1="contr.sum")}. See the \code{contrasts.arg} of \code{\link{model.matrix.default}}.}

\item{family}{family object}

\item{se}{(logical) compute standard error?}

\item{call}{original \code{glmmTMB} call}

\item{verbose}{whether progress indication should be printed to the console.}

\item{ziPredictCode}{zero-inflation code}

\item{doPredict}{flag to enable sds of predictions}

\item{whichPredict}{which observations in model frame represent predictions}

\item{REML}{whether to use REML estimation rather than maximum likelihood.}

\item{start}{starting values, expressed as a list with possible components \code{beta}, \code{betazi}, \code{betadisp} (fixed-effect parameters for conditional, zero-inflation, dispersion models); \code{b}, \code{bzi}, \code{bdisp} (conditional modes for conditional, zero-inflation, and dispersion models); \code{theta}, \code{thetazi}, \code{thetadisp} (random-effect parameters, on the standard deviation/Cholesky scale, for conditional, z-i, and disp models); \code{psi} (extra family parameters, e.g., shape for Tweedie models).}

\item{map}{a list specifying which parameter values should be fixed to a constant value rather than estimated. \code{map} should be a named list containing factors corresponding to a subset of the internal parameter names (see \code{start} parameter). Distinct factor values are fitted as separate parameter values, \code{NA} values are held fixed: e.g., \code{map=list(beta=factor(c(1,2,3,NA)))} would fit the first three fixed-effect parameters of the conditional model and fix the fourth parameter to its starting value. In general, users will probably want to use \code{start} to specify non-default starting values for fixed parameters. See \code{\link[TMB]{MakeADFun}} for more details.}

\item{sparseX}{see \code{\link{glmmTMB}}}

\item{control}{control parameters, see \code{\link{glmmTMBControl}}.}

\item{old_smooths}{(optional) smooth components from a previous fit: used when constructing a new model structure for prediction
from an existing model. A list of smooths for each model component; each smooth has sm and re elements}

\item{priors}{see \code{\link{priors}}}
}
\description{
Extract info from formulas, reTrms, etc., format for TMB
}
\keyword{internal}
