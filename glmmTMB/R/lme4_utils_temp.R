#### COPIED from lme4, until next version of lme4 (>= 1.18-1.9000) gets to CRAN

asDf0 <- function(x,nx,id=FALSE) {
    xt <- x[[nx]]
    ss <- utils::stack(xt)
    ss$ind <- factor(as.character(ss$ind), levels = colnames(xt))
    ss$.nn <- rep.int(stats::reorder(factor(rownames(xt)), xt[[1]],
                              FUN = mean,sort = sort), ncol(xt))
    ## allow 'postVar' *or* 'condVar' names
    pv <- attr(xt,"postVar")
    if (is.null(pv)) {
        pv <- attr(xt,"condVar")
    }
    if (!is.null(pv)) {
        tmpfun <- function(pvi) {
            unlist(lapply(1:nrow(pvi), function(i) sqrt(pvi[i, i, ])))
        }
        if (!is.list(pv)) {
            ss$se <- tmpfun(pv)
        } else {
            ## rely on ordering when unpacking!
            ss$se <- unlist(lapply(pv,tmpfun))
        }
    }
    if (id) ss$id <- nx
    return(ss)
}

