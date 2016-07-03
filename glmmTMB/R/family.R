##' Extended family functions for glmmTMB
##' 
##' @aliases family_glmmTMB
##' @param link (character) link function
##' @return returns a list with (at least) components
##' \item{family}{length-1 character vector giving the family name}
##' \item{link}{length-1 character vector specifying the link function}
##' \item{variance}{a function of either 1 (mean) or 2 (mean and dispersion
##' parameter) arguments giving the predicted variance}
##' 
##' @export
nbinom2 <- function(link="log") {
    return(list(family="nbinom2",link=link,
           variance=function(mu,disp) {
               mu*(1+mu/disp)
           }))
}

#' @rdname nbinom2
#' @export
nbinom1 <- function(link="log") {
    return(list(family="nbinom1",link=link,
           variance=function(mu,disp) {
               mu*disp
           }))
}

## FIXME: add truncated families (use unconditional variance?),
## check for zeros on initialization

## similar to mgcv::betar(), but simplified (variance has two parameters
##  rather than retrieving a variable from the environment); initialize()
##  tests for legal response values
betar <- function(link="logit") {
    return(list(family="betar",link=link,
                variance=function(mu,phi) {
                    mu*(1-mu)/(1+phi)
                },
                initialize=expression({
                    if (any(y <= 0 | y >= 1)) 
                        stop("y values must be 0 < y < 1")
                })))
}

## fixme: better name?

#' List model options that glmmTMB knows about
#'
#' @note these are all the options that are \emph{defined} internally; they have not necessarily all been \emph{implemented} (FIXME!)
#' @param what (character) which type of model structure to report on
#' ("all","family","link","covstruct")
#' @export
getCapabilities <- function(what="all") {
    switch(what,
           all=list(family=.valid_family,link=.valid_link,
                    covstruct=.valid_covstruct),
           family=.valid_family,
           link=.valid_link,
           covstruct=.valid_covstruct,
           stop(sprintf("unknown option %s",what)))
}
