# calculate the rank of a structure
get_rank <- function(reStruc){
  if(reStruc[["blockCode"]]==9){
    p <- reStruc$blockSize
    nt <- reStruc$blockNumTheta
    rank <- (2*p + 1 - sqrt((2*p+1)^2 - 8*nt))/2
  } else
    rank <- 0
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
##' @return A named list with three components (\code{cond}, \code{zi}, \code{disp}). 
##'   Each component contains a list with one element per random effect term. For rr 
##'   structures, each element contains:
##' \itemize{
##'   \item \code{lv} matrix of latent variables with columns named \code{lv1}, \code{lv2},..
##'   \item \code{fl} matrix of factor loadings
##' }  
##'   Non-rr terms return \code{NULL}.
##' @export rrTerms
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
