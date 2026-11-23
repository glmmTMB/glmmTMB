## rank of a reStrucTrm block:
## 0 for non-rr structures; 
## for rr  blockSize p and blockNumTheta nt 
get_rank <- function(reStrucTrm) {
  if (reStrucTrm[["blockCode"]] != .valid_covstruct[["rr"]]) return(0)
  p  <- reStrucTrm$blockSize
  nt <- reStrucTrm$blockNumTheta
  rank <- (2*p + 1 - sqrt((2*p+1)^2 - 8*nt))/2
  return(rank)
}

# Extract rr components 
extract_rr_info <- function(reTrms, reStruc, par_nm, pl, parRep){
  
  cnms <- reTrms$cnms   ## list of (named) terms and X columns
  flist <- reTrms$flist ## list of grouping variables
  flist_asgn <- attr(flist, "assign")
  levs <- lapply(seq_along(flist_asgn), function(i) levels(flist[[flist_asgn[i]]]))
  nc <- vapply(reStruc, function(x) x$blockSize, numeric(1)) ## number of RE params per block
  
  #function to split b by the random effect terms
  split.bseq <- function(reStruc){
    nc <- vapply(reStruc, function(x) x$blockSize, numeric(1)) ## number of RE params per block
    nb <- vapply(reStruc, function(x) x$blockReps, numeric(1)) ## number of blocks per RE (may != nlevs in some cases)
    ### splitting the b's into their respective random effects
    nbseq <- rep.int(seq_along(nb), nb * nc)       ## splitting vector
    return(nbseq)
  }
  nbseq <- split.bseq(reStruc)      ## splitting vector
  ml.b <- split(pl[[paste0("b", par_nm)]], nbseq)
  ml <- ml.b
  
  for (i in seq_along(ml.b)) {
    ml[[i]] <- matrix(ml.b[[i]], ncol = nc[i], byrow = TRUE,
                      dimnames = list(levs[[i]], cnms[[i]]))
  }
  
  rank <- vapply(reStruc, get_rank, FUN.VALUE=numeric(1))
  rr_pars <- vector(mode="list", length(rank))
  names(rr_pars) <- names(rank)
  
  for(j in seq_along(rr_pars)){
    nlv <- rank[[j]]
    if(nlv>0){
      rr_pars[[j]]$lv <- as.matrix(ml[[j]][,1:nlv])
      colnames(rr_pars[[j]]$lv) <- paste0("lv", 1:nlv)
      rr_pars[[j]]$fl <- parRep[[paste0("fact_load", par_nm)]][[j]]
      rownames(rr_pars[[j]]$fl) <- cnms[[j]]
      colnames(rr_pars[[j]]$fl) <- paste0("fl", 1:nlv)
    }
  }
  return(rr_pars)
}

##' Extract  reduced-rank (rr) factor loadings and latent variables
##' 
##' Extracts the factor loadings and latent variables from reduced-rank 
##' random effect structures in a fitted \code{glmmTMB} model.
##'
##' @name rrTerms
##' @param x a fitted \code{glmmTMB} model
##' @return A named list with three components:
##' \item{cond}{a list with one element per random-effect term in the
##'   conditional model.}
##' \item{zi}{a list with one element per random-effect term in the
##'   zero-inflation model.}
##' \item{disp}{a list with one element per random-effect term in the
##'   dispersion model.}
##' For an \code{rr} term, the element is a list containing:
##' \itemize{
##'   \item \code{lv}: a matrix of latent variables, with one row per level of
##'     the grouping factor and columns named \code{lv1}, \code{lv2}, \dots
##'   \item \code{fl}: a matrix of factor loadings
##' }
##' Terms that are not \code{rr} terms  return \code{NULL}.
##' @export
##' @examples
##' \dontrun{
##'   # Fit model with reduced-rank structure
##'   m1_rr <- glmmTMB(abund ~ Species + rr(Species + 0|id, d = 1), data = spider_long)
##'   m1_rr_trms <- rrTerms(m1_rr)
##'   m1_rr_trms$cond[[1]]$lv   # latent variables
##'   m1_rr_trms$cond[[1]]$fl   # factor loadings
##' }
rrTerms <- function(x){
  
  reT <- x$modelInfo$reTrms
  reS <- x$modelInfo$reStruc

  parRep <- x$obj$env$report(x$fit$parfull) #factor loading matrix
  pl <- x$obj$env$parList(x$fit$par, x$fit$parfull) #b latent variables
  
  comp_lst <- c("cond", "zi", "disp")
  par_nms <- c("", "zi", "disp")
  rrTrms <- vector(mode="list", length(comp_lst))
  names(rrTrms) <- comp_lst
  for (i in seq_along(comp_lst)) {
    rrTrms[[i]] <- extract_rr_info(reT[[comp_lst[i]]], reS[[paste0(comp_lst[i], "ReStruc")]], 
                              par_nms[i], pl, parRep)
  }
  
  return(rrTrms)
  
}

##' Plot latent variables from rr term
##'
##' Plots latent variables from a reduced-rank random-effect term of a
##' fitted \code{glmmTMB} model. When \code{biplot = TRUE} the factor loadings
##' are added as labelled points (or arrows).
##'
##' @name rrPlot
##' @param object a fitted \code{glmmTMB} model 
##' @param which_rr a list that selects which rr term to plot, with elements
##'   \itemize{
##'     \item \code{component}: \code{"cond"}, \code{"zi"} or \code{"disp"}
##'     \item \code{n}: the nth rr term within that component (out of all rr terms)
##'     \item \code{lvs}: the latent variables to plot
##'   }
##' @param biplot if \code{TRUE}, add the factor loadings as labelled points
##' @param arrows if \code{TRUE} (and \code{biplot = TRUE}), draw the loadings as
##'   labelled arrows from the origin instead of labelled points
##' @param labels if \code{TRUE}, label the latent-variable points with the
##'   levels of the grouping factor
##' @param scale a numeric to rescale the factor loadings relative to the latent variables when
##'   \code{biplot = TRUE}
##' @param \dots further arguments passed to \code{\link{plot}}
##' @return a list of the plotted latent variables (\code{lv}),
##'   the selection used (\code{which_rr}), and, for a biplot, the loadings
##'   (\code{fl}) and the arrow scaling (\code{scale})
##' @seealso \code{\link{rrTerms}}
##' @importFrom graphics plot points abline arrows text
##' @export
##' @examples
##' \dontrun{
##'   # Fit model with reduced-rank structure
##'   m1_rr <- glmmTMB(abund ~ Species + rr(Species + 0|id, d = 1), data = spider_long)
##'   rrPlot(m1_rr)
##' }
rrPlot <- function(object,
                   which_rr = list(component = "cond", n = 1, lvs = c(1, 2)),
                   biplot = FALSE,
                   arrows = FALSE,
                   labels = FALSE,
                   scale = NULL,
                   main = NULL,
                   ...) {

  if (!is.list(which_rr)) {
    stop("'which_rr' must be a list, e.g. list(component = \"cond\", n = 1, lvs = c(1, 2))",
         call. = FALSE)
  }
  defaults <- list(component = "cond", n = 1, lvs = c(1, 2))
  sel <- defaults
  sel[names(which_rr)] <- which_rr
  component <- match.arg(sel$component, c("cond", "zi", "disp"))
  n   <- sel$n
  which_lvs <- sel$lvs
  if (length(n) != 1 || !is.numeric(n) || n < 1 || n != round(n)) {
    stop("'n' in 'which_rr' must be a positive integer", call. = FALSE)
  }

  ## get rr terms in th component (non-rr entries are NULL)
  all_rr_terms <- rrTerms(object)
  rr_comp <- all_rr_terms[[component]]
  rr <- rr_comp[!vapply(rr_comp, is.null, NA)]
  if (length(rr) == 0) {
    stop("no reduced-rank (rr) terms found in the '", component, "' component",
         call. = FALSE)
  }
  if (n > length(rr)) {
    stop("'", component, "' has ", length(rr), " rr term(s): ",
         paste0(seq_along(rr), " = ", names(rr), collapse = "; "),
         call. = FALSE)
  }

  lv_fl_n <- rr[[n]]
  ## check the latent variables to plot
  d <- ncol(lv_fl_n$lv)
  if (d < 2) {
    stop("this rr term has a single latent variable (d = 1): nothing to plot against",
         call. = FALSE)
  }
  if (!is.numeric(which_lvs) || any(which_lvs < 1) || any(which_lvs > d)) {
    stop("'lvs' in 'which_rr' must be between 1 and ", d, call. = FALSE)
  }

  lv <- lv_fl_n$lv[, which_lvs, drop = FALSE]
  fl <- lv_fl_n$fl[, which_lvs, drop = FALSE]
  k <- if(biplot) { if(is.null(scale)) max(abs(lv)) / max(abs(fl)) else scale } else { 1 }

  ## common limits so that both axes have the same scale and loadings fit
  lim <- range(lv, if (biplot) fl * k)
  lim <- lim + c(-1, 1) * 0.1 * diff(lim)

  plot(lv, xlim = lim, ylim = lim,  pch = 16, ...)
  if (labels) text(lv, labels = rownames(lv), pos = 3, cex = 0.7)

  if (biplot) {
    if (arrows) {
      arrows(0, 0, fl[, 1] * k, fl[, 2] * k, length = 0.08, col = "blue")
    } 
    text(fl[, 1] * k, fl[, 2] * k, labels = rownames(fl),
         pos = 3, col = "blue", cex = 0.8)
  }
  invisible(c(list(lv = lv,
                   which_rr = list(component = component, n = n, lvs = which_lvs)),
              if (biplot) list(fl = fl, scale = k)))
}