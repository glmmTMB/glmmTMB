#' Compute likelihood profiles for a fitted model
#'
#' @inheritParams confint.glmmTMB
#' @param fitted a fitted \code{glmmTMB} object
#' @param parm which parameters to profile, specified
#' \itemize{
#' \item by index (position)
#' \item by name (matching the row/column names of \code{vcov(object,full=TRUE)})
#' \item as \code{"theta_"} (random-effects variance-covariance parameters) or \code{"beta_"} (conditional and zero-inflation parameters)
#' }
#' @param level_max maximum confidence interval target for profile
#' @param npts target number of points in (each half of) the profile (\emph{approximate})
#' @param stepfac initial step factor (fraction of estimated standard deviation)
#' @param trace print tracing information? If \code{trace=FALSE} or 0,
#' no tracing; if \code{trace=1}, print names of parameters currently
#' being profiled; if \code{trace>1}, turn on tracing for the
#' underlying \code{\link{tmbprofile}} function
#' @param stderr standard errors to use as a scaling factor when picking step
#' sizes to compute the profile; by default (if \code{stderr} is
#' \code{NULL}, or \code{NA} for a particular element),
#' uses the estimated (Wald) standard errors of the parameters
#' @param ... additional arguments passed to \code{\link{tmbprofile}}
#' @return An object of class \code{profile.glmmTMB}, which is also a
#' data frame, with columns \code{.par} (parameter being profiled),
#' \code{.focal} (value of focal parameter), value (negative log-likelihood).
#' @importFrom stats profile
#' @examples
#' \dontrun{
#' m1 <- glmmTMB(count~ mined + (1|site),
#'        zi=~mined, family=poisson, data=Salamanders)
#' salamander_prof1 <- profile(m1, parallel="multicore",
#'                             ncpus=2, trace=1)
#' ## testing
#' salamander_prof1 <- profile(m1, trace=1,parm=1)
#' salamander_prof1M <- profile(m1, trace=1,parm=1, npts = 4)
#' salamander_prof2 <- profile(m1, parm="theta_")
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
                            stderr = NULL,
                            trace = FALSE,
                            parallel = c("no", "multicore", "snow"),
                            ncpus = getOption("profile.ncpus", 1L),
                            cl = NULL,
                            ...) {

    if (isREML(fitted)) stop("can't compute profiles for REML models at the moment (sorry)")
    plist <- parallel_default(parallel,ncpus)
    parallel <- plist$parallel
    do_parallel <- plist$do_parallel

    trace <- as.numeric(trace)

    ytol <- qchisq(level_max,1)
    ystep <- ytol/npts

    ## don't suppress sigma profiling (full=TRUE)
    parm <- getParms(parm, fitted, full=TRUE)

    ## only need selected SDs
    sds <- sqrt(diag(vcov(fitted,full=TRUE)))
    sds <- sds[parm]

    if (!is.null(stderr)) {
        if (length(stderr) != length(sds)) {
            if (length(stderr)==1) {
                sds <- rep(stderr,length(sds))
            } else {
                stop(
          sprintf("length(stderr) should equal 1 or number of parameters (%d)",
                  length(parm)))
            }
        } else {
            sds[!is.na(stderr)] <- stderr[!is.na(stderr)]
        }
    }

    if (any(sds>1e3)) {
        warning("very large standard errors for parameters: ",
                names(sds)[sds>1e3])
    }

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


    FUN <- local({
        function(p,s) {
            if (trace>0) cat("parameter",p,"\n")
            n_orig <- openmp(NULL)
            openmp(n = fitted$modelInfo$parallel)
            on.exit(openmp(n_orig))
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
    ## stringsAsFactors=TRUE below to maintain backward compatibility/
    ## ordering of rows in confint.profile.glmmTMB ...
    dfun <- function(x,n) {
        dd0 <- data.frame(n,x, stringsAsFactors=TRUE)
        names(dd0)[1:2] <- c(".par",".focal")
        dd0$value <- dd0$value - min(dd0$value,na.rm=TRUE)
        return(dd0)
    }
    dd <- Map(dfun, L, names(sds))
    dd <- do.call(rbind,c(dd, list(stringsAsFactors=TRUE)))
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
    qval <- 0.5*qchisq(level,df=1)
    ci_fun <- function(dd) {
        dd <- dd[!duplicated(dd$.focal),] ## unique values: WHY??
        hf <- with(dd,factor(.focal>.focal[which.min(value)],
                   levels=c("FALSE","TRUE")))
        halves <- split(dd,hf)
        res <- vapply(halves,ci_fun_half,numeric(1))
        a <- (1 - level)/2
        a <- c(a, 1 - a)
        names(res) <- format.perc(a, 3)
        return(res)
    }
    ## fit spline and invert for one half (lower, upper) of the profile
    ci_fun_half <- function(hh) {
        if (nrow(hh)==0) return(NA_real_)
        if (max(hh$value,na.rm=TRUE)<qval) {
            restr_prof_flag <- TRUE
        }
        for_spl <- splines::interpSpline(value~.focal,hh)
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
