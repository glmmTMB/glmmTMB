\name{Owls}
\title{Begging by Owl Nestlings}
\docType{data}
\alias{Owls}
\alias{OwlModel}
\alias{OwlModel_nb1_bs}
\alias{OwlModel_nb1_bs_mcmc}
\description{
  Begging by owl nestlings
}
\usage{data(Owls)}
\format{
  The \code{Owls} data set is a data frame with
  599 observations on the following variables:
  \describe{
    \item{\code{Nest}}{a factor describing individual nest locations}
    \item{\code{FoodTreatment}}{(factor) food treatment: \code{Deprived}
      or \code{Satiated}}
    \item{\code{SexParent}}{(factor) sex of provisioning parent: \code{Female} or \code{Male}}
    \item{\code{ArrivalTime}}{a numeric vector}
    \item{\code{SiblingNegotiation}}{a numeric vector}
    \item{\code{BroodSize}}{brood size}
    \item{\code{NegPerChick}}{number of negotations per chick}
  }
}
\source{
  Roulin, A. and L. Bersier (2007)
  Nestling barn owls beg more intensely in the presence of their mother
  than in the presence of their father.
  \emph{Animal Behaviour} \bold{74} 1099--1106.
  \doi{10.1016/j.anbehav.2007.01.027};
  \url{http://www.highstat.com/Books/Book2/ZuurDataMixedModelling.zip}
}
\references{
  Zuur, A. F., E. N. Ieno, N. J. Walker, A. A. Saveliev, and G. M. Smith (2009)
  \emph{Mixed Effects Models and Extensions in Ecology with R}; Springer.
}
\note{Access to data kindly provided by Alain Zuur}
\examples{
data(Owls, package = "glmmTMB")
require("lattice")
bwplot(reorder(Nest,NegPerChick) ~ NegPerChick | FoodTreatment:SexParent,
       data=Owls)
dotplot(reorder(Nest,NegPerChick) ~ NegPerChick| FoodTreatment:SexParent,
        data=Owls)
\dontrun{
## Fit negative binomial model with "constant" Zero Inflation :
owls_nb1 <- glmmTMB(SiblingNegotiation ~ FoodTreatment*SexParent +
                                    (1|Nest)+offset(log(BroodSize)),
              family = nbinom1(), zi = ~1, data=Owls)
owls_nb1_bs <- update(owls_nb1,
                      . ~ . - offset(log(BroodSize)) + log(BroodSize))
fixef(owls_nb1_bs)
}
}
\keyword{datasets}
