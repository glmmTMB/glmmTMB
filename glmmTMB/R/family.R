##' Family functions for glmmTMB
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


##' @export
genpois <- function(link="log") {
    return(list(family="genpois",link=link,
           variance=function(mu,disp) {
               mu*disp
           }))
}

##' @export
compois <- function(link="log") {
    return(list(family="compois",link=link,
           variance=function(mu,disp) {
               stop("variance for compois family not yet implemented")
               #mu*disp #FIXME: need to call C++ calc_loglambda
           }))
}

#' List model options that glmmTMB knows about
#'
#' @note these are all the options that are \emph{defined} internally; they have not necessarily all been \emph{implemented} (FIXME!)
#' @param what (character) which type of model structure to report on
#' ("all","family","link","covstruct")
#' @param check (logical) do brute-force checking to test whether families are really implemented (only available for \code{what="family"})
#' @return if \code{check==FALSE}, returns a vector of the names (or a list of name vectors) of allowable entries; if \code{check==TRUE}, returns a logical vector of working families
#' 
#' @export
getCapabilities <- function(what="all",check=FALSE) {
    if (!check) {
        switch(what,
               all=lapply(list(family=.valid_family,link=.valid_link,
                        covstruct=.valid_covstruct),names),
               family=names(.valid_family),
               link=names(.valid_link),
               covstruct=names(.valid_covstruct),
               stop(sprintf("unknown option %s",what)))
    } else {
        ## run dummy models to see if we get a family-not-implemented error
        if (what!="family") stop("'check' option only available for families")
        families <- names(.valid_family)
        family_OK <- setNames(rep(TRUE,length(.valid_family)),families)
        y <- 1:3 ## dummy
        for (f in families) {
            tt1 <- utils::capture.output(tt0 <- suppressMessages(suppressWarnings(
                try(glmmTMB(y~1,
                            family=list(family=f,link="identity")),
                    silent=TRUE))))
            family_OK[f] <- !(inherits(tt0,"try-error") &&
                              grepl("Family not implemented!",tt0))

        }
        return(family_OK)
    }
}
