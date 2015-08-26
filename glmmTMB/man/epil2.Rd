\name{epil2}
\alias{epil2}
\docType{data}
\title{Seizure Counts for Epileptics}
\description{
  Extended version of the \code{epil} dataset of the \pkg{MASS} package.
  The three transformed variables \code{Visit}, \code{Base}, and
  \code{Age} used by Booth et al. (2003) have been added to \code{epil}.
}
\usage{epil2}
\format{
  A data frame with 236 observations on the following 12 variables:
  \describe{
    \item{\code{y}}{an integer vector.}
    \item{\code{trt}}{a factor with levels \code{"placebo"} and
      \code{"progabide"}.}
    \item{\code{base}}{an integer vector.}
    \item{\code{age}}{an integer vector.}
    \item{\code{V4}}{an integer vector.}
    \item{\code{subject}}{an integer vector.}
    \item{\code{period}}{an integer vector.}
    \item{\code{lbase}}{a numeric vector.}
    \item{\code{lage}}{a numeric vector.}
    \item{Visit}{\code{(rep(1:4,59) - 2.5) / 5}.}
    \item{Base}{\code{log(base/4)}.}
    \item{Age}{\code{log(age)}.}
  }
}
\references{
  Booth, J.G., G. Casella, H. Friedl, and J.P. Hobert. 2003. Negative
  binomial loglinear mixed models. Statistical Modelling 3:179--191.
}
\seealso{
  \code{\link{glmmadmb}}.

  \code{\link{glmmADMB-package}}.
}
\examples{
epil2$subject <- factor(epil2$subject)
\dontrun{
(fm <- glmmadmb(y~Base*trt+Age+Visit+(Visit|subject),
                data=epil2, family="nbinom"))
}
}
\keyword{datasets}
