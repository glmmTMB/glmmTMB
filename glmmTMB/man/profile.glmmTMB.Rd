% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/profile.R
\name{profile.glmmTMB}
\alias{profile.glmmTMB}
\alias{confint.profile.glmmTMB}
\title{Compute likelihood profiles for a fitted model}
\usage{
\method{profile}{glmmTMB}(
  fitted,
  parm = NULL,
  level_max = 0.99,
  npts = 8,
  stepfac = 1/4,
  stderr = NULL,
  trace = FALSE,
  parallel = c("no", "multicore", "snow"),
  ncpus = getOption("profile.ncpus", 1L),
  cl = NULL,
  ...
)

\method{confint}{profile.glmmTMB}(object, parm = NULL, level = 0.95, ...)
}
\arguments{
\item{fitted}{a fitted \code{glmmTMB} object}

\item{parm}{which parameters to profile, specified
\itemize{
\item by index (position)
\item by name (matching the row/column names of \code{vcov(object,full=TRUE)})
\item as \code{"theta_"} (random-effects variance-covariance parameters) or \code{"beta_"} (conditional and zero-inflation parameters)
}}

\item{level_max}{maximum confidence interval target for profile}

\item{npts}{target number of points in (each half of) the profile (\emph{approximate})}

\item{stepfac}{initial step factor (fraction of estimated standard deviation)}

\item{stderr}{standard errors to use as a scaling factor when picking step
sizes to compute the profile; by default (if \code{stderr} is
\code{NULL}, or \code{NA} for a particular element),
uses the estimated (Wald) standard errors of the parameters}

\item{trace}{print tracing information? If \code{trace=FALSE} or 0,
no tracing; if \code{trace=1}, print names of parameters currently
being profiled; if \code{trace>1}, turn on tracing for the
underlying \code{\link[TMB]{tmbprofile}} function}

\item{parallel}{method (if any) for parallel computation}

\item{ncpus}{number of CPUs/cores to use for parallel computation}

\item{cl}{cluster to use for parallel computation}

\item{...}{additional arguments passed to \code{\link[TMB]{tmbprofile}}}

\item{object}{a fitted profile (\code{profile.glmmTMB}) object}

\item{level}{confidence level}
}
\value{
An object of class \code{profile.glmmTMB}, which is also a
data frame, with columns \code{.par} (parameter being profiled),
\code{.focal} (value of focal parameter), value (negative log-likelihood).
}
\description{
Compute likelihood profiles for a fitted model
}
\details{
Fits natural splines separately to the points from each half of the profile for each
specified parameter (i.e., values above and below the MLE), then finds the inverse functions
to estimate the endpoints of the confidence interval
}
\examples{
\dontrun{
m1 <- glmmTMB(count~ mined + (1|site),
       zi=~mined, family=poisson, data=Salamanders)
salamander_prof1 <- profile(m1, parallel="multicore",
                            ncpus=2, trace=1)
## testing
salamander_prof1 <- profile(m1, trace=1,parm=1)
salamander_prof1M <- profile(m1, trace=1,parm=1, npts = 4)
salamander_prof2 <- profile(m1, parm="theta_")

}
salamander_prof1 <- readRDS(system.file("example_files","salamander_prof1.rds",package="glmmTMB"))
if (require("ggplot2")) {
    ggplot(salamander_prof1,aes(.focal,sqrt(value))) +
        geom_point() + geom_line()+
        facet_wrap(~.par,scale="free_x")+
    geom_hline(yintercept=1.96,linetype=2)
}
salamander_prof1 <- readRDS(system.file("example_files","salamander_prof1.rds",package="glmmTMB"))
confint(salamander_prof1)
confint(salamander_prof1,level=0.99)
}
