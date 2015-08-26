\name{Owls}
\alias{Owls}
\alias{OwlModel}
\alias{OwlModel_nb1_bs}
\alias{OwlModel_nb1_bs_mcmc}

\docType{data}
\title{
Begging by owl nestlings
}
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
  For convenience in testing methods, several pre-fitted models are provided:
  \itemize{
    \item \code{OwlModel} is the results of fitting the model
      \code{glmmadmb(SiblingNegotiation~FoodTreatment*SexParent+
	(1|Nest)+offset(log(BroodSize)),family="nbinom",
	zeroInflation=TRUE,data=Owls)}
    \item \code{OwlModel_nb1_bs}: \code{glmmadmb(SiblingNegotiation~FoodTreatment*SexParent+BroodSize+
      (1|Nest),family="nbinom1",zeroInflation=TRUE,data=Owls)}
    \item \code{OwlModel_nb1_bs_mcmc}: \code{glmmadmb(SiblingNegotiation~FoodTreatment*SexParent+BroodSize+
      (1|Nest),family="nbinom1",zeroInflation=TRUE,MCMC=TRUE,data=Owls)}
    }
}
\source{
  Roulin, A. and L. Bersier. 2007. Nestling barn owls beg more intensely
  in the presence of their mother than in the presence of their father. Animal
  Behaviour 74:1099-1106. \url{http://www.sciencedirect.com/science/article/B6W9W-4PK8B6H-8/2/e43cfbaad4dc0bb2207adfc54a460c89};
  \url{http://www.highstat.com/Book2/ZuurDataMixedModelling.zip}
}
\references{
Zuur, A. F., E. N. Ieno, N. J. Walker, A. A. Saveliev, and
G. M. Smith. 2009. Mixed Effects Models and Extensions in Ecology with R. Springer.
}
\note{Access to data kindly provided by Alain Zuur}
\examples{
data(Owls)
library(lattice)
bwplot(reorder(Nest,NegPerChick)~NegPerChick|FoodTreatment:SexParent,
data=Owls)
dotplot(reorder(Nest,NegPerChick)~NegPerChick|FoodTreatment:SexParent,
data=Owls)
getdata("OwlModel")
summary(OwlModel)
}
\keyword{datasets}
