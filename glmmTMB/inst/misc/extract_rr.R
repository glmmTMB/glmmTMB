extract_rr <- function(object){
  listname <- "cond"
  cnms <- object$modelInfo$reTrms[[listname]]$cnms   ## list of (named) terms and X columns
  flist <- object$modelInfo$reTrms[[listname]]$flist ## list of grouping variables
  flist_asgn <- attr(flist, "assign")
  levs <- lapply(seq_along(flist_asgn), function(i) levels(flist[[flist_asgn[i]]]))
  
  reStruc <- object$modelInfo$reStruc[[paste0(listname, "ReStruc")]] ## random-effects structure
  nc <- vapply(reStruc, function(x) x$blockSize, numeric(1)) ## number of RE params per block
  pl <- object$obj$env$parList(object$fit$par, object$fit$parfull)
  
  #function to split b by the random effect terms
  split.bseq <- function(object){
    listname <- "cond"
    reStruc <- object$modelInfo$reStruc[[paste0(listname, "ReStruc")]] ## random-effects structure
    nc <- vapply(reStruc, function(x) x$blockSize, numeric(1)) ## number of RE params per block
    nb <- vapply(reStruc, function(x) x$blockReps, numeric(1)) ## number of blocks per RE (may != nlevs in some cases)
    ### splitting the b's into their respective random effects
    nbseq <- rep.int(seq_along(nb), nb * nc)       ## splitting vector
    return(nbseq)
  }
  nbseq <- split.bseq(object)      ## splitting vector
  ml.b <- split(pl$b, nbseq)
  ml <- ml.b
  
  for (i in seq_along(ml.b)) {
    ml[[i]] <- matrix(ml.b[[i]], ncol = nc[i], byrow = TRUE,
                      dimnames = list(levs[[i]], cnms[[i]]))
  }
  
  get_rank <- function(x){
    if(x[["blockCode"]]==9){
      p <- x$blockSize
      nt <- x$blockNumTheta
      rank <- (2*p + 1 - sqrt((2*p+1)^2 - 8*nt))/2
    } else
      rank <- 0
    return(rank)
  }
  
  rank <- vapply(object$modelInfo$reStruc$condReStruc,
                 get_rank,
                 FUN.VALUE=numeric(1))
  nlv <- rank[rank > 0]
  rrName <- names(nlv)
  rrBlock <- which(rank > 0)
  b = ml[[rrBlock]][,1:nlv]
  colnames(b) <- paste0("lv", 1:nlv)
  fact_load <- object$obj$env$report(object$fit$parfull)$fact_load[[rrBlock]]
  rownames(fact_load) <- cnms[[rrBlock]]
  
  return(list(fl = fact_load, b = b))
}
