#' Compute likelihood profiles for a fitted model
#' 
#' @param fitted a fitted \code{glmmTMB} object
#' @param which which parameters to profile, specified by position
#' @param parallel method (if any) for parallel computation
#' @param ncpus number of CPUs/cores to use for parallel computation
#' @param cl cluster to use for parallel computation
#' @param trace print tracing information? If \code{trace=FALSE} or 0,
#' no tracing; if \code{trace=1}, print names of parameters currently
#' being profiled; if \code{trace>1}, turn on tracing for the
#' underlying \code{\link{tmbprofile}} function
#' @param stepfac initial step factor (fraction of estimated standard deviation)
#' @param ... additional arguments passed to \code{\link{tmbprofile}}
#' @return An object of class \code{profile.glmmTMB}, which is also a
#' data frame, with columns \code{.par} (parameter being profiled),
#' \code{.focal} (value of focal parameter), value (negative log-likelihood).
#' @examples
#' m1 <- glmmTMB(count~ mined + (1|site), 
#'        zi=~mined, family=poisson, data=Salamanders)
#' \dontrun{
#' salamander_prof1 <- profile(m1,ytol=0.5, parallel="multicore",
#'                             ncpus=2, trace=1)
#' }
#' salamander_prof1 <- readRDS(system.file("example_files","salamander_prof1.rds",package="glmmTMB"))
#' if (require("ggplot2")) {
#'     ggplot(salamander_prof1,aes(.focal,value)) + geom_point() + geom_line()+
#'         facet_wrap(~.par,scale="free_x")
#' }
#' @importFrom TMB tmbprofile
#' @export
profile.glmmTMB <- function(fitted, which=NULL,
                            parallel = c("no", "multicore", "snow"),
                            ncpus = getOption("profile.ncpus", 1L),
                            cl = NULL,
                            trace = FALSE,
                            stepfac = 1/4,
                            ...) {

    ## lots of boilerplate parallel-handling stuff, copied from lme4
    if (missing(parallel)) parallel <- getOption("profile.parallel", "no")
    parallel <- match.arg(parallel)
    trace <- as.numeric(trace)
    do_parallel <- (parallel != "no" && ncpus > 1L)
    if (do_parallel && parallel == "multicore" &&
        .Platform$OS.type == "windows") {
        warning("no multicore on Windows, falling back to non-parallel")
        parallel <- "no"
    }

    ## get sds
    vv <- vcov(fitted,full=TRUE)
    sds <- sqrt(diag(vv))

    ## get pars: for now, all of them
    pars <- seq_along(sds)

    stds <- sqrt
    FUN <- local({
        function(p,s) {
            if (trace>0) cat("parameter",p,"\n")
            return(tmbprofile(fitted$obj,name=p,h=s/4,
                              trace=(trace>1),...))
        }
    })
    if (do_parallel) {
        if (parallel == "multicore") {
            L <- parallel::mcmapply(FUN, pars, sds, mc.cores = ncpus,
                                    SIMPLIFY=FALSE)
        } else if (parallel=="snow") {
            if (is.null(cl)) {
                ## start cluster
                new_cl <- TRUE
                cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
            }
            ## run
            L <- parallel::clusterMap(cl, FUN, pars, sds)
            if (new_cl) {
                ## stop cluster
                parallel::stopCluster(cl)
            }
        } else {
            L <- mapply(FUN, pars, sds)
        }
    }
    dfun <- function(x,n) {
        dd0 <- data.frame(.par=n,x)
        names(dd0)[2] <- ".focal"
        return(dd0)
    }
    dd <- Map(dfun, L, colnames(vv)[pars])
    dd <- do.call(rbind,dd)
    class(dd) <- c("profile.glmmTMB","data.frame")
    return(dd)
}

#' @importFrom splines interpSpline backSpline
confint.profile.glmmTMB <- function(object, parm=NULL, level = 0.95, ...) {
    ## FIXME: lots of bulletproofing:
    ##   non-monotonic values: error and/or linear interpolation
    ##   non-monotonic spline,
    ## find CIs for a single parameter
    ci_fun <- function(dd) {
        dd <- dd[!duplicated(dd$.focal),] ## unique values: WHY??
        dd$min <- min(dd$value)
        halves <- with(dd,split(dd,.focal>.focal[which.min(value)]))
        res <- vapply(halves,ci_fun_half,numeric(1))
        names(res) <- c("lwr","upr")
        return(res)
    }
    ## fit spline and invert for one half (lower, upper) of the profile
    ci_fun_half <- function(hh) {
        for_spl <- splines::interpSpline(value-min~.focal,hh)
        bak_spl <- splines::backSpline(for_spl)
        predict(bak_spl,c(1.96))$y
    }
    objList <- split(object,object$.par)
    if (is.null(parm)) {
        parm <- seq_along(objList)
    } else {
        stop("parameter selection not yet implemented")
    }
    ci_mat <- t(vapply(objList,ci_fun,numeric(2)))
    return(ci_mat)
}
