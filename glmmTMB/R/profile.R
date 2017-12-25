#' Compute likelihood profiles for a fitted model
#' 
#' @param fitted a fitted \code{glmmTMB} object
#' @param parm which parameters to profile, specified by index (position)
#' @param level_max maximum confidence interval target for profile
#' @param npts target number of points in (each half of) the profile (\emph{approximate})
#' @param stepfac initial step factor (fraction of estimated standard deviation)
#' @param trace print tracing information? If \code{trace=FALSE} or 0,
#' no tracing; if \code{trace=1}, print names of parameters currently
#' being profiled; if \code{trace>1}, turn on tracing for the
#' underlying \code{\link{tmbprofile}} function
#' @param parallel method (if any) for parallel computation
#' @param ncpus number of CPUs/cores to use for parallel computation
#' @param cl cluster to use for parallel computation
#' @param ... additional arguments passed to \code{\link{tmbprofile}}
#' @return An object of class \code{profile.glmmTMB}, which is also a
#' data frame, with columns \code{.par} (parameter being profiled),
#' \code{.focal} (value of focal parameter), value (negative log-likelihood).
#' @examples
#' \dontrun{
#' m1 <- glmmTMB(count~ mined + (1|site), 
#'        zi=~mined, family=poisson, data=Salamanders)
#' salamander_prof1 <- profile(m1, parallel="multicore",
#'                             ncpus=2, trace=1)
#' ## testing
#' salamander_prof1 <- profile(m1, trace=1,parm=1)
#' salamander_prof1M <- profile(m1, trace=1,parm=1, npts = 4)
#' 
#' }
#' salamander_prof1 <- readRDS(system.file("example_files","salamander_prof1.rds",package="glmmTMB"))
#' if (require("ggplot2")) {
#'     ggplot(salamander_prof1,aes(.focal,sqrt(value))) +
#'         geom_point() + geom_line()+
#'         facet_wrap(~.par,scale="free_x")+
#'     geom_hline(yintercept=1.96,linetype=2)
#' }
#' @importFrom TMB tmbprofile
#' @export
profile.glmmTMB <- function(fitted,
                            parm=NULL,
                            level_max = 0.99,
                            npts = 8,
                            stepfac = 1/4,
                            trace = FALSE,
                            parallel = c("no", "multicore", "snow"),
                            ncpus = getOption("profile.ncpus", 1L),
                            cl = NULL,
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

    ytol <- qchisq(level_max,1)
    ystep <- ytol/npts

    vv <- vcov(fitted,full=TRUE)
    sds <- sqrt(diag(vv))

    ## get pars: need to match up names with internal positions
    if (is.null(parm)) parm <-  seq_along(sds)
    if (FALSE) {
        ## would like complete solution for assigning names to components
        ## (cond (fix/theta), zi (fix/theta), disp (fix/theta))
        ## and matching order of parameters in object ...
        nn <- names(object$obj$par)
        drop_null <- function(x) x[lengths(x)>0]
        ff <- drop_null(fixef(object))
        pn <- lapply(ff,names)
        nm <- unlist(mapply(paste,names(ff),pn,MoreArgs=list(sep="_"),
                  SIMPLIFY=FALSE))
    }

    ## only need selected SDs
    sds <- sds[parm]

    FUN <- local({
        function(p,s) {
            if (trace>0) cat("parameter",p,"\n")
            return(tmbprofile(fitted$obj,
                              name=p,
                              h=s/4,
                              ytol=ytol,
                              ystep=ystep,
                              trace=(trace>1),...))
        }
    })
    if (do_parallel) {
        if (parallel == "multicore") {
            L <- parallel::mcmapply(FUN, parm, sds, mc.cores = ncpus,
                                    SIMPLIFY=FALSE)
        } else if (parallel=="snow") {
            if (is.null(cl)) {
                ## start cluster
                new_cl <- TRUE
                cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
            }
            ## run
            L <- parallel::clusterMap(cl, FUN, parm, sds)
            if (new_cl) {
                ## stop cluster
                parallel::stopCluster(cl)
            }
        }
    } else { ## non-parallel
        L <- Map(FUN, parm, sds)
    }
    dfun <- function(x,n) {
        dd0 <- data.frame(n,x)
        names(dd0)[1:2] <- c(".par",".focal")
        dd0$value <- dd0$value - min(dd0$value,na.rm=TRUE)
        return(dd0)
    }
    dd <- Map(dfun, L, colnames(vv)[parm])
    dd <- do.call(rbind,dd)
    class(dd) <- c("profile.glmmTMB","data.frame")
    return(dd)
}

#' @param object a fitted profile (\code{profile.glmmTMB}) object
#' @param level confidence level
#' @rdname profile.glmmTMB
#' @details Fits natural splines separately to the points from each half of the profile for each
#' specified parameter (i.e., values above and below the MLE), then finds the inverse functions
#' to estimate the endpoints of the confidence interval
#' @examples
#' salamander_prof1 <- readRDS(system.file("example_files","salamander_prof1.rds",package="glmmTMB"))
#' confint(salamander_prof1)
#' confint(salamander_prof1,level=0.99)
#' @importFrom splines interpSpline backSpline
#' @export
confint.profile.glmmTMB <- function(object, parm=NULL, level = 0.95, ...) {
    ## FIXME: lots of bulletproofing:
    ##   non-monotonic values: error and/or linear interpolation
    ##   non-monotonic spline,
    ## find CIs for a single parameter

    qval <- qnorm((1+level)/2)
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
        if (max(sqrt(hh$value),na.rm=TRUE)<qval) {
            restr_prof_flag <- TRUE
        }
        for_spl <- splines::interpSpline(sqrt(value)~.focal,hh)
        bak_spl <- splines::backSpline(for_spl)
        predict(bak_spl,qval)$y
    }
    objList <- split(object,object$.par)
    if (is.null(parm)) {
        parm <- seq_along(objList)
    }
    restr_prof_flag <- FALSE
    ci_mat <- t(vapply(objList[parm],ci_fun,numeric(2)))
    if (restr_prof_flag) {
        warning("max profile values less than target: consider increasing level_max when computing profiles")
    }
    return(ci_mat)
}
