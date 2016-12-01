\name{Salamanders}
\title{Repeated counts of salamanders in streams}
\alias{Salamanders}
\docType{data}
\description{
  A dataset containing counts of salamanders with site covariates and sampling covariates. Each of 23 sites were sampled 4 times.
  When using this data, please cite Price et al. (2016) as well as the the Dryad data package (Price et al. 2015). 
}
\usage{data(Salamanders)}
\format{
  A data frame with 644 observations on the following 10 variables:
  \describe{
	\item{site}{name of a location where repeated samples were taken}
	\item{mined}{factor indicating whether the site was affected by mountian top removal coal mining} 
	\item{cover}{amount of cover objects in the stream (scaled)} 
	\item{sample}{repeated sample}
	\item{DOP}{Days since precipitation (scaled)}
	\item{Wtemp}{water temperature (scaled)}
	\item{DOY}{day of year (scaled)}
	\item{spp}{abreviated species name, possibly also life stage}
	\item{count}{number of observed salamanders}  
}
}
\references{
    Price SJ, Muncy BL, Bonner SJ, Drayer AN, Barton CD (2016) Effects of mountaintop removal mining and valley filling on the occupancy and abundance of stream salamanders. \emph{Journal of Applied Ecology} \bold{53} 459--468. \url{http://dx.doi.org/10.1111/1365-2664.12585}
  
    Price SJ, Muncy BL, Bonner SJ, Drayer AN, Barton CD (2015) Data from: Effects of mountaintop removal mining and valley filling on the occupancy and abundance of stream salamanders. \emph{Dryad Digital Repository}. \url{http://dx.doi.org/10.5061/dryad.5m8f6} 
}

\examples{
require("glmmTMB")
data(Salamanders)
zipm3 = glmmTMB(count~spp * mined + (1|site), zi=~spp * mined, Salamanders, family="poisson")
}

\keyword{datasets}
