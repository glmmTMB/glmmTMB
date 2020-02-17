# added 2017-12-13 by J. Fox
# 2017-12-14: improved recovery of model data
#             removed faulty one-step approximations
# 2018-01-28: fix computation of Cook's D for lme models
# 2018-05-23: fixed bug when more than one grouping variable (reported by Maarten Jung)
# 2018-06-07: skip plot of "sigma^2" in GLMM if dispersion fixed to 1; improved labelling for covariance components
# 2018-11-04: tweak to dfbetas.influence.merMod() suggested by Ben Bolker.
# 2018-11-09: parallel version of influence.merMod()
# 2020-02-17: works with tibble data (Ben Bolker/Joseph O'Brien)

## influence diagnostics for mixed models

## copied from car (unexported, don't want to rely on car:::combineLists)
combineLists <- function(..., fmatrix="list", flist="c", fvector="rbind", 
                         fdf="rbind", recurse=FALSE){
    # combine lists of the same structure elementwise
    
    # ...: a list of lists, or several lists, each of the same structure
    # fmatrix: name of function to apply to matrix elements
    # flist: name of function to apply to list elements
    # fvector: name of function to apply to data frame elements
    # recurse: process list element recursively
    
    frecurse <- function(...){
        combineLists(..., fmatrix=fmatrix, fvector=fvector, fdf=fdf, 
                     recurse=TRUE)
    }
    
    if (recurse) flist="frecurse"
    list.of.lists <- list(...)
    if (length(list.of.lists) == 1){
        list.of.lists <- list.of.lists[[1]]
        list.of.lists[c("fmatrix", "flist", "fvector", "fdf")] <- 
            c(fmatrix, flist, fvector, fdf)
        return(do.call("combineLists", list.of.lists))
    }
    if (any(!sapply(list.of.lists, is.list))) 
        stop("arguments are not all lists")
    len <- sapply(list.of.lists, length)
    if (any(len[1] != len)) stop("lists are not all of the same length")
    nms <- lapply(list.of.lists, names)
    if (any(unlist(lapply(nms, "!=", nms[[1]])))) 
        stop("lists do not all have elements of the same names")
    nms <- nms[[1]]
    result <- vector(len[1], mode="list")
    names(result) <- nms
    for(element in nms){
        element.list <- lapply(list.of.lists, "[[", element)
#        clss <- sapply(element.list, class)
        clss <- lapply(element.list, class)
#        if (any(clss[1] != clss)) stop("list elements named '", element,
        if (!all(vapply(clss, function(e) all(e == clss[[1L]]), NA)))
          stop("list elements named '", element, "' are not all of the same class")
        
        is.df <- is.data.frame(element.list[[1]])
        fn <- if (is.matrix(element.list[[1]])) fmatrix 
        else if (is.list(element.list[[1]]) && !is.df) flist 
        else if (is.vector(element.list[[1]])) fvector
        else if (is.df) fdf
        else stop("list elements named '", element, 
                  "' are not matrices, lists, vectors, or data frames")
        result[[element]] <- do.call(fn, element.list)
    }
    result
}

## copied from lme4 
namedList <- function (...)  {
    L <- list(...)
    snm <- sapply(substitute(list(...)), deparse)[-1]
    if (is.null(nm <- names(L))) 
        nm <- snm
    if (any(nonames <- nm == "")) 
        nm[nonames] <- snm[nonames]
    setNames(L, nm)
}

##' compute influence measures for [g]lmerMod (from lme4) or glmmTMB (from glmmTMB) fitted models
##' @param model fitted model
##' @param groups (character) names of group(s) to leave out/compute sensitivity for: ".case"  means observation-level sensitivity
##' @param data data frame to use in refitting (will be taken from the model call if possible, but it may be necessary to specify this if the influence measures are being computed in a new environment)
##' @param maxfun maximum number of function evaluations
##' @param ncores number of parallel cores to use for computation
##' @param component for glmmTMB models, which component to use (conditional, zero-inflated, or dispersion)
##' @examples
##' source(system.file("other_methods","influence_mixed.R",package="glmmTMB"))
##' m1 <- glmmTMB(count ~ mined + (1|site),
##'   zi=~mined,
##'   family=poisson, data=Salamanders)
##' i1 <- influence_mixed(m1,groups="site")
##' ## check that it works with tibbles too ...
##' m2 <- update(m1, data=tibble::as_tibble(Salamanders))
##' i2 <- influence_mixed(m2,groups="site")
##' if (require(car)) {
##'    stopifnot(all.equal(dfbeta(i1),dfbeta(i2)))
##' }
influence_mixed <- function(model, groups=".case", data, maxfun=1000,
                            ncores=getOption("mc.cores", 2L),
                            component=c("cond","zi","disp"),
                            progress=FALSE,
                            ...) {
    component <- match.arg(component)
    if (inherits(model,"glmmTMB")) {
        fe <- function(x) fixef(x)[[component]]
        VC <- function(x) VarCorr(x)[[component]]
        .vcov <- function(x) Matrix::as.matrix(vcov(x)[[component]])
    } else {
        fe <- fixef
        VC <- VarCorr
        .vcov <- function(x) Matrix::as.matrix(vcov(x))
    }
    
    if (is.infinite(ncores)) {
        ncores <- parallel::detectCores(logical=FALSE)
    }
    if (missing(data)){
        data <- getCall(model)$data
        data <- if (!is.null(data)) eval(data, parent.frame())
        else stop("model did not use the data argument")
    }
    if (".case" %in% groups) {
        data$.case <- rownames(data)
    }
    else if (length(groups) > 1){
        del.var <- paste0(groups, collapse=".")
        ## Reduce(interaction,groups) ?
        data[[del.var]] <- apply(data, 1, paste0, collapse=".")
        groups <- del.var
    }
    unique.del <- unique(data[[groups]])
    data[[".groups"]] <- data[[groups]]
    par <- list(theta=getME(model, "theta"))
    if (inherits(model, "glmerMod") || inherits(model,"glmmTMB")) {
        par$beta <- fe(model)
    }
    fixed <- fe(model)
    Vs <- VC(model)
    nms <- names(Vs)
    sep <- ":"
    if (length(nms) == 1) {
        nms <- ""
        sep <- ""
    }
    vc <- sigma(model)^2
    names(vc) <- "sigma^2"
    for (i in 1:length(Vs)){
        V <- Vs[[i]]
        c.names <- colnames(V)
        e.names <- outer(c.names, c.names, function(a, b) paste0("C[", a, ",", b, "]"))
        diag(e.names) <- paste0("v[", c.names, "]")
        v <- V[lower.tri(V, diag=TRUE)]
        names(v) <- paste0(nms[i], sep, e.names[lower.tri(e.names, diag=TRUE)])
        vc <- c(vc, v)
    }
    if (inherits(model,"lmerMod")) {
        control <- lme4::lmerControl(optCtrl=list(maxfun=maxfun))
    } else if (inherits(model, "glmerMod")) {
        control <- lme4:: glmerControl(optCtrl=list(maxfun=maxfun))
    } else if (inherits(model, "glmmTMB")) {
        ## ???
        control <- glmmTMBControl(optCtrl=list(iter.max=maxfun,eval.max=maxfun))
    }
    deleteGroup <- function(del){
        if (progress) cat(".")
        data_sub <- data[data$.groups!=del,]
        mod.1 <- suppressWarnings(update(model, data=data_sub,
                                         start=par,
                                         control=control))
        ## platform-independent?
        if (inherits(mod.1,"glmmTMB")) {
            feval <- mod.1$fit$evaluations[["function"]]
            converged <- mod.1$fit$convergence==0
        } else { ## [g]lmerMod
            opt <- mod.1@optinfo
            feval <- opt$feval
            converged <- opt$conv$opt == 0 && length(opt$warnings) == 0
        }
        fixed.1 <- fe(mod.1)
        Vs.1 <- VC(mod.1)
        vc.0 <- sigma(mod.1)^2
        for (V in Vs.1){
            vc.0 <- c(vc.0, V[lower.tri(V, diag=TRUE)])
        }
        vc.1 <- vc.0
        vcov.1 <<- .vcov(mod.1)
        namedList(fixed.1, vc.1, vcov.1, converged, feval)
    }
    result <- if(ncores >= 2){
        message("using a cluster of ", ncores, " cores")
        cl <- parallel::makeCluster(ncores)
        on.exit(parallel::stopCluster(cl))
        parallel::clusterExport(cl, c("namedList", "combineLists"))
        parallel::clusterEvalQ(cl, require("lme4"))
        parallel::clusterEvalQ(cl, require("glmmTMB"))
        parallel::clusterApply(cl, unique.del, deleteGroup)
    } else {
        lapply(unique.del, deleteGroup)
    }
    result <- combineLists(result)
    fixed.1 <- result$fixed.1
    rownames(fixed.1) <- unique.del
    colnames(fixed.1) <- names(fixed)
    vc.1 <- result$vc.1
    rownames(vc.1) <- unique.del
    colnames(vc.1) <- names(vc)
    feval <- as.vector(result$feval)
    converged <- as.vector(result$converged)
    vcov.1 <- result$vcov.1
    names(vcov.1) <- names(feval) <- names(converged) <- unique.del
    left <- "[-"
    right <- "]"
    if (groups == ".case") {
        groups <- "case"
    }
    nms <- c("fixed.effects", paste0("fixed.effects", left, groups, right),
             "var.cov.comps", paste0("var.cov.comps", left, groups, right),
             "vcov", paste0("vcov", left, groups, right),
             "groups", "deleted", "converged", "function.evals")
    result <- list(fixed, fixed.1, vc, vc.1,
                   .vcov(model), vcov.1, groups, unique.del, converged, feval)
    names(result) <- nms
    ## call this influence.merMod since all methods are currently written for it
    class(result) <- "influence.merMod"
    result
}
