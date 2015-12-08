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


