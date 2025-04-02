data("HawCon", package = "qra")
library(glmmTMB)
library(splines)
## devtools::load_all("~/R/pkgs/glmmTMB/glmmTMB")

## load/process data from qra package
CCnum <- match("CommonName", names(HawCon))
names(HawCon)[CCnum] <- "CN"
HawCon <- within(HawCon, {
  trtGp <- factor(paste0(CN,LifestageTrt, sep=":"))
  trtGpRep <- paste0(CN,LifestageTrt,":",RepNumber)
  scTime <- scale(TrtTime)
  obs <- factor(1:nrow(HawCon))
})

form <- cbind(Dead,Live)~0+trtGp/TrtTime+(1|trtGpRep)
HCbb.cll <- glmmTMB::glmmTMB(form, dispformula=~trtGp+ns(scTime,2), 
                    family=glmmTMB::betabinomial(link="cloglog"),
                    data=HawCon)
vv0 <- vcov(HCbb.cll)
HCbin.cll <- update(HCbb.cll, family=binomial(link="cloglog"))
vv1 <- vcov(HCbin.cll)
