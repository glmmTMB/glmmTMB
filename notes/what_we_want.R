library(glmmTMB)
data(Owls, package="glmmADMB")
Owls <- transform(Owls, ArrivalTime=c(scale(ArrivalTime,scale=FALSE)),
                  NCalls=SiblingNegotiation)
tmbz <- glmmTMB(NCalls ~ FoodTreatment + ArrivalTime + (ArrivalTime|Nest) + (1|SexParent),
                ziformula = ~1 + (1|Nest) + (1|SexParent),
                data=Owls, family=poisson(link="log"))
object <- tmbz

################################################################################

pl <- object$obj$env$parList(object$fit$par, object$obj$env$last.par.best)
Z <- getME(object, "Z")
recond <- structure(pl$b, names=colnames(Z))
Zzi <- getME(object, "Zzi")
rezi <- structure(pl$bzi, names=colnames(Zzi))

grpVar <- object$modelInfo$grpVar

output <- list(conditional_model=list(), zero_inflation=list())

## conditional_model
cnms <- object$modelInfo$reTrms$condList$cnms
flist <- object$modelInfo$reTrms$condList$flist

## cond[1]
nterms <- length(cnms[[1]])
nlev <- nlevels(flist[[1]])
beg <- 1
end <- beg + nterms * nlev - 1

m <- matrix(recond[beg:end], ncol=nterms, byrow=TRUE)
rownames(m) <- levels(flist[[1]])
colnames(m) <- cnms[[1]]
output$conditional_model[[grpVar[1]]] <- m

## cond[2]
nterms <- length(cnms[[2]])
nlev <- nlevels(flist[[2]])
beg <- end + 1
end <- beg + nterms * nlev - 1

m <- matrix(recond[beg:end], ncol=nterms, byrow=TRUE)
rownames(m) <- levels(flist[[2]])
colnames(m) <- cnms[[2]]
output$conditional_model[[grpVar[2]]] <- m

## zero_inflation
cnms <- object$modelInfo$reTrms$ziList$cnms
flist <- object$modelInfo$reTrms$ziList$flist

## zero[1]
nterms <- length(cnms[[1]])
nlev <- nlevels(flist[[1]])
beg <- 1
end <- beg + nterms * nlev - 1

m <- matrix(rezi[beg:end], ncol=nterms, byrow=TRUE)
rownames(m) <- levels(flist[[1]])
colnames(m) <- cnms[[1]]
output$zero_inflation[[grpVar[1]]] <- m

## zero[2]

nterms <- length(cnms[[2]])
nlev <- nlevels(flist[[2]])
beg <- end + 1
end <- beg + nterms * nlev - 1

m <- matrix(rezi[beg:end], ncol=nterms, byrow=TRUE)
rownames(m) <- levels(flist[[2]])
colnames(m) <- cnms[[2]]
output$zero_inflation[[grpVar[2]]] <- m

## done
