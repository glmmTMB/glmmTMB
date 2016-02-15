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
