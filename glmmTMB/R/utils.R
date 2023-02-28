## backward compat (copied from lme4)
if((Rv <- getRversion()) < "3.2.1") {
    lengths <- function (x, use.names = TRUE) vapply(x, length, 1L, USE.NAMES = use.names)
}
rm(Rv)

## generate a list with names equal to values
## See also: \code{tibble::lst}, \code{Hmisc::llist}
namedList <- function (...) {
    L <- list(...)
    snm <- sapply(substitute(list(...)), deparse)[-1]
    if (is.null(nm <- names(L)))
        nm <- snm
    if (any(nonames <- nm == ""))
        nm[nonames] <- snm[nonames]
    setNames(L, nm)
}

## Sparse Schur complement (Marginal of precision matrix)
##' @importFrom Matrix Cholesky solve
GMRFmarginal <- function(Q, i, ...) {
    ind <- seq_len(nrow(Q))
    i1 <- (ind)[i]
    i0 <- setdiff(ind, i1)
    if (length(i0) == 0)
        return(Q)
    Q0 <- as(Q[i0, i0, drop = FALSE], "symmetricMatrix")
    L0 <- Cholesky(Q0, ...)
    ans <- Q[i1, i1, drop = FALSE] - Q[i1, i0, drop = FALSE] %*%
        solve(Q0, as.matrix(Q[i0, i1, drop = FALSE]))
    ans
}

parallel_default <- function(parallel=c("no","multicore","snow"),ncpus=1) {
    ##  boilerplate parallel-handling stuff, copied from lme4
    if (missing(parallel)) parallel <- getOption("profile.parallel", "no")
    parallel <- match.arg(parallel)
    do_parallel <- (parallel != "no" && ncpus > 1L)
    if (do_parallel && parallel == "multicore" &&
        .Platform$OS.type == "windows") {
        warning("no multicore on Windows, falling back to non-parallel")
        parallel <- "no"
    }
    return(list(parallel=parallel,do_parallel=do_parallel))
}

##' translate vector of correlation parameters to correlation values
##' @param theta vector of internal correlation parameters (elements of scaled Cholesky factor, in \emph{row-major} order)
##' @return a vector of correlation values (\code{get_cor}) or glmmTMB scaled-correlation parameters (\code{put_cor})
##' @details These functions follow the definition at \url{http://kaskr.github.io/adcomp/classdensity_1_1UNSTRUCTURED__CORR__t.html}:
##' if \eqn{L} is the lower-triangular matrix with 1 on the diagonal and the correlation parameters in the lower triangle, then the correlation matrix is defined as \eqn{\Sigma = D^{-1/2} L L^\top D^{-1/2}}{Sigma = sqrt(D) L L' sqrt(D)}, where \eqn{D = \textrm{diag}(L L^\top)}{D = diag(L L')}. For a single correlation parameter \eqn{\theta_0}{theta0}, this works out to \eqn{\rho = \theta_0/\sqrt{1+\theta_0^2}}{rho = theta0/sqrt(1+theta0^2)}. The \code{get_cor} function returns the elements of the lower triangle of the correlation matrix, in column-major order.
##' @examples
##' th0 <- 0.5
##' stopifnot(all.equal(get_cor(th0),th0/sqrt(1+th0^2)))
##' get_cor(c(0.5,0.2,0.5))
##' C <- matrix(c(1,  0.2,  0.1,
##'              0.2,  1, -0.2,
##'              0.1,-0.2,   1),
##'            3, 3)
##' ## test: round-trip (almostl results in lower triangle only)
##' stopifnot(all.equal(get_cor(put_cor(C)),
##'                    C[lower.tri(C)]))
##' @export
get_cor <- function(theta) {
  n <- as.integer(round(0.5 * (1 + sqrt(1 + 8 * length(theta)))))
  R <- diag(n)
  R[upper.tri(R)] <- theta
  R[] <- crossprod(R) # R <- t(R) %*% R
  scale <- 1 / sqrt(diag(R))
  R[] <- scale * R * rep(scale, each = n) # R <- cov2cor(R)
  R[lower.tri(R)]
}

##' @rdname get_cor
##' @param C a correlation matrix
##' @export
put_cor <- function(C) {
    cc <- chol(C)
    cc2 <- t(cc %*% diag(1/diag(cc)))
    cc2[lower.tri(cc2)]
}

hasRandom <- function(x) {
    pl <- getParList(x)
    return(length(unlist(pl[grep("^theta",names(pl))]))>0)
}

##' retrieve parameters by name or index
##' @param parm parameter specifier
##' @param object fitted glmmTMB object
##' @param full include simple dispersion parameter?
##' @param include_mapped include mapped parameter indices?
##' @noRd
getParms <- function(parm=NULL, object, full=FALSE, include_mapped = FALSE) {
    vv <- vcov(object, full=TRUE, include_mapped = include_mapped)
    sds <- sqrt(diag(vv))
    pnames <- names(sds) <- rownames(vv)       ## parameter names (user-facing)
    ee <- object$obj$env

    ## don't use object$obj$env$random; we want to keep "beta" vals, which may be
    ## counted as "random" if using REML
    drop_rand <- function(x) x[!x %in% c("b", "bzi")]
    if (!include_mapped) {
        intnames <- drop_rand(names(ee$last.par))
    } else {
        pl <- ee$parList()
        intnames <- drop_rand(rep(names(pl), lengths(pl)))
    }
    if (length(pnames) != length(sds)) { ## shouldn't happen ...
        stop("length mismatch between internal and external parameter names")
    }

    if (is.null(parm)) {
        if (!full && trivialDisp(object)) {
            parm <- grep("betad", intnames, invert=TRUE)
        } else {
            parm <- seq_along(sds)
        }
    }
    if (is.character(parm)) {
        if (identical(parm,"psi_")) {
            parm <- grep("^psi",intnames)
        } else if (identical(parm,"theta_")) {
            parm <- grep("^theta",intnames)
        } else if (identical(parm,"beta_")) {
            if (trivialDisp(object)) {
                ## include conditional and zi params
                ##   but not dispersion params
                parm <- grep("^beta(zi)?$",intnames)
            } else {
                parm <- grep("beta",intnames)
            }
        } else if (identical(parm, "disp_") ||
                   identical(parm, "sigma")) {
            parm <- grep("^betad", intnames)
        } else { ## generic parameter vector
            nparm <- match(parm,pnames)
            if (any(is.na(nparm))) {
                stop("unrecognized parameter names: ",
                     parm[is.na(nparm)])
            }
            parm <- nparm
        }
    }
    return(parm)
}

isREML <- function(x) {
    if (is.null(REML <- x$modelInfo$REML)) {
        ## let vcov work with old (pre-REML option) stored objects
        REML <- FALSE
    }
    return(REML)
}

## action: message, warning, stop
check_dots <- function(..., .ignore = NULL, .action="stop") {
    L <- list(...)
    if (length(.ignore)>0) {
        L <- L[!names(L) %in% .ignore]
    }
    if (length(L)>0) {
        FUN <- get(.action)
        FUN("unknown arguments: ",
            paste(names(L), collapse=","))
    }
    return(NULL)
}

if (getRversion()<"4.0.0") {
    deparse1 <- function (expr, collapse = " ", width.cutoff = 500L, ...) {
        paste(deparse(expr, width.cutoff, ...), collapse = collapse)
    }
}

## in case these are useful, we can document and export them later ...
#' @importFrom stats rnbinom qnbinom dnbinom pnbinom

rnbinom1 <- function(n, mu, phi) {
    ## var = mu*(1+phi) = mu*(1+mu/k) -> k = mu/phi
    rnbinom(n, mu=mu, size=mu/phi)
}

dnbinom1 <- function(x, mu, phi, ...) {
    dnbinom(x, mu=mu, size=mu/phi, ...)
}

pnbinom1 <- function(q, mu, phi, ...) {
    pnbinom(q, mu=mu, size=mu/phi, ...)
}

qnbinom1 <- function(p, mu, phi, ...) {
    qnbinom(p, mu=mu, size=mu/phi, ...)
}

nullSparseMatrix <- function() {
    argList <- list(
        dims=c(0,0),
        i=integer(0),
        j=integer(0),
        x=numeric(0))
    if (utils::packageVersion("Matrix")<"1.3.0") {
        do.call(Matrix::sparseMatrix, c(argList, list(giveCsparse=FALSE)))
    } else {
        do.call(Matrix::sparseMatrix, c(argList, list(repr="T")))
    }
}

#' Check for version mismatch in dependent binary packages
#' @param dep_pkg upstream package
#' @param this_pkg downstream package
#' @param write_file (logical) write version file and quit?
#' @param warn give warning?
#' @return logical: TRUE if the binary versions match
#' @importFrom utils packageVersion
#' @export
checkDepPackageVersion <- function(dep_pkg = "TMB",
                                   this_pkg = "glmmTMB",
                                   write_file = FALSE,
                                   warn = TRUE) {
    cur_dep_version <- as.character(packageVersion(dep_pkg))
    fn <- sprintf("%s-version", dep_pkg)
    if (write_file) {
        cat(sprintf("current %s version=%s: writing file\n", dep_pkg, cur_dep_version))
        writeLines(cur_dep_version, con = fn)
        return(cur_dep_version)
    }
    fn <- system.file(fn, package=this_pkg)
    built_dep_version <- scan(file=fn, what=character(), quiet=TRUE)
    result_ok <- identical(built_dep_version, cur_dep_version)
    if(warn && !result_ok) {
        warning(
            "Package version inconsistency detected.\n",
            sprintf("%s was built with %s version %s",
                    this_pkg, dep_pkg, built_dep_version),
            "\n",
            sprintf("Current %s version is %s",
                    dep_pkg, cur_dep_version),
            "\n",
            sprintf("Please re-install %s from source ", this_pkg),
            "or restore original ",
            sQuote(dep_pkg), " package (see '?reinstalling' for more information)"
        )
    }
    return(result_ok)
}

#' @name reinstalling
#' @rdname reinstalling
#' @title Reinstalling binary dependencies
#'
#' @description The \code{glmmTMB} package depends on several upstream packages, which it
#' uses in a way that depends heavily on their internal (binary) structure.
#' Sometimes, therefore, installing an update to one of these packages will
#' require that you re-install a \emph{binary-compatible} version of \code{glmmTMB},
#' i.e. a version that has been compiled with the updated version of the upstream
#' package.
#' \itemize{
#' \item If you have development tools (compilers etc.) installed, you
#' should be able to re-install a binary-compatible version of the package by running
#' \code{install.packages("glmmTMB", type="source")}. If you want to install
#' the development version of \code{glmmTMB} instead, you can use
#' \code{remotes::install_github("glmmTMB/glmmTMB/glmmTMB")}.
#' (On Windows, you can install development tools following the instructions at
#' \url{https://cran.r-project.org/bin/windows/Rtools/}; on MacOS, see
#' \url{https://mac.r-project.org/tools/}.)
#'
#' \item If you do \emph{not} have development tools and can't/don't want to
#' install them (and so can't install packages with compiled code from source),
#' you have two choices:
#' \itemize{
#' \item revert the upstream package(s) to their previous binary version. For example, using the
#' \code{checkpoint} package:
#' \preformatted{
#' ## load (installing if necessary) the checkpoint package
#' while (!require("checkpoint")) install.packages("checkpoint")
#' ## retrieve build date of installed version of glmmTMB
#' bd <- as.character(asDateBuilt(
#'       packageDescription("glmmTMB",fields="Built")))
#' oldrepo <- getOption("repos")
#' use_mran_snapshot(bd) ## was setSnapshot() pre-checkpoint v1.0.0
#' install.packages("TMB")
#' options(repos=oldrepo) ## restore original repo
#' }
#' A similar recipe (substituting \code{Matrix} for \code{TMB} and \code{TMB} for \code{glmmTMB})
#' can be used if you get warnings about an incompatibility between \code{TMB} and \code{Matrix}.
#' \item hope that the glmmTMB maintainers have posted a binary
#' version of the package that works with your system; try installing it via
#' \code{install.packages("glmmTMB",repos="https://glmmTMB.github.io/glmmTMB/repos",type="binary")}
#' If this doesn't work, please file an issue (with full details about your
#' operating system and R version) asking the maintainers to build and
#' post an appropriate binary version of the package.
#' }
#' }
NULL

#' Check OpenMP status
##'
##' Checks whether OpenMP has been successfully enabled for this
##' installation of the package. (Use the \code{parallel} argument
##' to \code{\link{glmmTMBControl}}, or set \code{options(glmmTMB.cores=[value])},
##' to specify that computations should be done in parallel.)
##' @seealso \code{\link[TMB]{benchmark}}, \code{\link{glmmTMBControl}}
##' @return \code{TRUE} or {FALSE} depending on availability of OpenMP
##' @export
omp_check <- function() {
    .Call("omp_check", PACKAGE="glmmTMB")
}

get_pars <- function(object, unlist=TRUE) {
    ee <- object$obj$env
    x <- ee$last.par.best
    ## work around built-in default to parList, which
    ##  is bad if no random component
    if (length(ee$random)>0) x <- x[-ee$random]
    p <- ee$parList(x=x)
    if (!unlist) return(p)
    p <- unlist(p[names(p)!="b"])  ## drop primary RE
    names(p) <- gsub("[0-9]+$","",names(p)) ## remove disambiguators
    return(p)
}


## replacement for (unexported) TMB:::isNullPointer
isNullPointer <- function(x) {
    attributes(x) <- NULL
    identical(x, new("externalptr"))
}

#' conditionally update glmmTMB object fitted with an old TMB version
#'
#' @rdname gt_load
#' @param oldfit a fitted glmmTMB object
#' @export
up2date <- function(oldfit) {
  openmp(1)  ## non-parallel/make sure NOT grabbing all the threads!
  if (isNullPointer(oldfit$obj$env$ADFun$ptr)) {
      obj <- oldfit$obj
      ee <- obj$env
      if ("thetaf" %in% names(ee$parameters)) {
          ee$parameters$psi <- ee$parameters$thetaf
          ee$parameters$thetaf <- NULL
          pars <- c(grep("last\\.par", names(ee), value = TRUE),
                    "par")
          for (p in pars) {
              if (!is.null(nm <- names(ee[[p]]))) {
                  names(ee[[p]])[nm == "thetaf"] <- "psi"
              }
          }
      }
      ee2 <- oldfit$sdr$env
      if ("thetaf" %in% names(ee2$parameters)) {
          ee2$parameters$psi <- ee2$parameters$thetaf
          ee2$parameters$thetaf <- NULL
      }
      oldfit$obj <- with(ee,
                       TMB::MakeADFun(data,
                                      parameters,
                                      map = map,
                                      random = random,
                                      silent = silent,
                                      DLL = "glmmTMB"))
      oldfit$obj$env$last.par.best <- ee$last.par.best
  }
  return(oldfit)
}



#' Load data from system file, updating glmmTMB objects
#'
#' @param fn partial path to system file (e.g. test_data/foo.rda)
#' @param verbose print names of updated objects?
#' @param mustWork fail if file not found?
#' @export
gt_load <- function(fn, verbose=FALSE, mustWork = FALSE) {
    sf <- system.file(fn, package = "glmmTMB")
    found_file <- file.exists(sf)
    if (mustWork && !found_file) {
        stop("couldn't find system file ", sf)
    }

    L <- load(sf)
    for (m in L) {
        if (inherits(get(m), "glmmTMB")) {
            if (verbose) cat(m,"\n")
            assign(m, up2date(get(m)))
        }
        assign(m, get(m), parent.env(), envir = parent.frame())
    }
    return(found_file)
}

#' truncated distributions
#'
#' Probability functions for k-truncated Poisson and negative binomial distributions.
#' @param x value
#' @param size number of trials/overdispersion parameter
#' @param mu mean parameter
#' @param k truncation parameter
#' @param log (logical) return log-probability?
#' @export
dtruncated_nbinom2 <- function(x, size, mu, k=0, log=FALSE) {
    y <- ifelse(x<=k,-Inf,
                dnbinom(x, mu=mu, size=size, log=TRUE) -
                pnbinom(k, mu=mu, size=size, lower.tail=FALSE,
                        log.p=TRUE))
    if (log) return(y) else return(exp(y))
}

#' @rdname dtruncated_nbinom2
#' @param lambda mean parameter
#' @importFrom stats dpois
#' @export
dtruncated_poisson <- function(x,lambda,k=0,log=FALSE) {
    y <- ifelse(x<=k,-Inf,
                dpois(x,lambda,log=TRUE) -
                ppois(k, lambda=lambda, lower.tail=FALSE,
                      log.p=TRUE))
    if (log) return(y) else return(exp(y))
}

#' @rdname dtruncated_nbinom2
#' @param phi overdispersion parameter
#' @export
dtruncated_nbinom1 <- function(x, phi, mu, k=0, log=FALSE) {
    ## V=mu*(1+phi) = mu*(1+mu/k) -> k=mu/phi
    size <- mu/phi
    y <- ifelse(x<=k,-Inf,
                dnbinom(x,mu=mu, size=size,log=TRUE) -
                pnbinom(k, mu=mu, size=size, lower.tail=FALSE,
                        log.p=TRUE))
    if (log) return(y) else return(exp(y))
}


## utilities for constructing lists of parameter names

## for matching map names vs nameList components ...
par_components <- c("beta","betazi","betad","theta","thetazi","psi")

getAllParnames <- function(object, full) {
                           
  mkNames <- function(tag) {
      X <- getME(object,paste0("X",tag))
      if (trivialFixef(nn <- colnames(X),tag)
          ## if 'full', keep disp even if trivial, if used by family
          && !(full && tag =="d" &&
               (usesDispersion(family(object)$family) && !zeroDisp(object)))) {
          return(character(0))
      }
      return(paste(tag,nn,sep="~"))
  }

  nameList <- setNames(list(colnames(getME(object,"X")),
                       mkNames("zi"),
                       mkNames("d")),
                names(cNames))

  if(full) {
      ## FIXME: haven't really decided if we should drop the
      ##   trivial variance-covariance dispersion parameter ??
      ## if (trivialDisp(object))
      ##    res <- covF[-nrow(covF),-nrow(covF)]

      reNames <- function(tag) {
        re <- object$modelInfo$reStruc[[paste0(tag,"ReStruc")]]
        num_theta <- vapply(re,"[[","blockNumTheta", FUN.VALUE = numeric(1))
        nn <- mapply(function(n,L) paste(n, seq(L), sep="."),
                     names(re), num_theta)
        if (length(nn) == 0) return(nn)
        return(paste("theta",gsub(" ", "", unlist(nn)), sep="_"))
      }
      ## nameList for estimated variables;
      nameList <- c(nameList,
                    list(theta = reNames("cond"), thetazi = reNames("zi")))

      ##
      if (length(fp <- family_params(object)) > 0) {
          nameList <- c(nameList, list(psi = names(fp)))
      }
      
  }

    return(nameList)
}


getEstParnames <- function(object, full) {
    nameList <- getAllParnames(object, full)
    map <- object$obj$env$map
    if (length(map)>0) {
        ## fullNameList for all variables, including mapped vars
        ## (nameList will get reduced shortly)
        for (m in seq_along(map)) {
            if (length(NAmap <- which(is.na(map[[m]])))>0) {
                w <- match(names(map)[m],par_components) ##
                if (length(nameList)>=w) { ## may not exist if !full
                    nameList[[w]] <- nameList[[w]][-NAmap]
                }
            }
        }
    }
    return(nameList)
}

## OBSOLETE: delete eventually

## reassign predvars to have term vars in the right order,
##  but with 'predvars' values inserted where appropriate
fix_predvars <- function(pv,tt) {
    if (length(tt)==3) {
        ## convert two-sided to one-sided formula
        tt <- RHSForm(tt, as.form=TRUE)
    }
    ## ugh, deparsing again ...
    tt_vars <- vapply(attr(tt, "variables"), deparse1, character(1))[-1]
    ## remove terminal paren - e.g. match term poly(x, 2) to
    ##   predvar poly(x, 2, <stuff>)
    ## beginning of string, including open-paren, colon
    ##  but not *first* comma nor arg ...
    ##  could possibly try init_regexp <- "^([^,]+).*" ?
    init_regexp <- "^([(^:_.[:alnum:]]+).*"
    tt_vars_short <- gsub(init_regexp,"\\1",tt_vars)
    if (is.null(pv) || length(tt_vars)==0) return(NULL)
    new_pv <- quote(list())
    ## maybe multiple variables per pv term ... [-1] ignores head
    ## FIXME: test for really long predvar strings ????
    pv_strings <- vapply(pv,deparse1,FUN.VALUE=character(1))[-1]
    pv_strings <- gsub(init_regexp,"\\1",pv_strings)
    for (i in seq_along(tt_vars)) {
        w <- match(tt_vars_short[[i]],pv_strings)
        if (!is.na(w)) {
            new_pv[[i+1]] <- pv[[w+1]]
        } else {
            ## insert symbol from term vars
            new_pv[[i+1]] <- as.symbol(tt_vars[[i]])
        }
    }
    return(new_pv)
}

make_pars <- function(pars, ...) {
    ## FIXME: check for name matches, length matches etc.
    L <- list(...)
    for (nm in names(L)) {
        pars[names(pars) == nm] <- L[[nm]]
    }
    return(pars)
}

##' Simulate from covariate/metadata in the absence of a real data set (EXPERIMENTAL)
##'
##' See \code{vignette("sim", package = "glmmTMB")} for more details and examples,
##' and \code{vignette("covstruct", package = "glmmTMB")}
##' for more information on the parameterization of different covariance structures.
##' 
##' @param object a \emph{one-sided} model formula (e.g. \code{~ a + b + c}
##' (peculiar naming is for consistency with the generic function, which typically
##' takes a fitted model object)
##' @param nsim number of simulations
##' @param seed random-number seed
##' @param newdata a data frame containing all variables listed in the formula,
##' \emph{including} the response variable (which needs to fall within
##' the domain of the conditional distribution, and should probably not
##' be all zeros, but whose value is otherwise irrelevant)
##' @param newparams a list of parameters containing sub-vectors
##' (\code{beta}, \code{betazi}, \code{betad}, \code{theta}, etc.) to
##' be used in the model
##' @param ... other arguments to \code{glmmTMB} (e.g. \code{family})
##' @param show_pars (logical) print structure of parameter vector and stop without simulating?
##' @examples
##' ## use Salamanders data for structure/covariates
##' simulate_new(~ mined + (1|site),
##'              zi = ~ mined,
##'              newdata = Salamanders, show_pars  = TRUE)
##' sim_count <- simulate_new(~ mined + (1|site),
##'              newdata = Salamanders,
##'              zi = ~ mined,
##'              family = nbinom2,
##'              newparams = list(beta = c(2, 1),
##'                          betazi = c(-0.5, 0.5), ## logit-linear model for zi
##'                          betad = log(2), ## log(NB dispersion)
##'                          theta = log(1)) ## log(among-site SD)
##' )
##' head(sim_count[[1]])
##' @export
simulate_new <- function(object,
                         nsim = 1,
                         seed = NULL,
                         newdata, newparams, ..., show_pars = FALSE) {
    if (!is.null(seed)) set.seed(seed)
    ## truncate
    if (length(object) == 3) stop("simulate_new should take a one-sided formula")
    ## fill in fake LHS
    form <- object
    form[[3]] <- form[[2]]
    form[[2]] <- quote(..y)
    ## insert a legal value: 1.0 is OK as long as family != "beta_family"
    newdata[["..y"]] <- if (!identical(list(...)$family, "beta_family")) 1.0 else 0.5
    r1 <- glmmTMB(form,
              data = newdata,
              ...,
              doFit = FALSE)
## construct TMB object, but don't fit it
    r2 <- fitTMB(r1, doOptim = FALSE)
    if (show_pars) return(r2$env$last.par)
    pars <- do.call("make_pars",
                    c(list(r2$env$last.par), newparams))
    replicate(nsim, r2$simulate(par = pars)$yobs, simplify = FALSE)
}
