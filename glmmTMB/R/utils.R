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

## from length of (strictly) lower triangle, compute dimension of matrix
get_matdim <- function(ntri) {
    as.integer(round(0.5 * (1 + sqrt(1 + 8 * ntri))))
}

##' translate vector of correlation parameters to correlation values
##' @param theta vector of internal correlation parameters (elements of scaled Cholesky factor, in \emph{row-major} order)
##' @param return_val return a vector of correlation values from the lower triangle ("vec"), or the full correlation matrix ("mat")? 
##' @return a vector of correlation values (\code{get_cor}) or glmmTMB scaled-correlation parameters (\code{put_cor})
##' @details These functions follow the definition at \url{http://kaskr.github.io/adcomp/classdensity_1_1UNSTRUCTURED__CORR__t.html}:
##' if \eqn{L} is the lower-triangular matrix with 1 on the diagonal and the correlation parameters in the lower triangle, then the correlation matrix is defined as \eqn{\Sigma = D^{-1/2} L L^\top D^{-1/2}}{Sigma = sqrt(D) L L' sqrt(D)}, where \eqn{D = \textrm{diag}(L L^\top)}{D = diag(L L')}. For a single correlation parameter \eqn{\theta_0}{theta0}, this works out to \eqn{\rho = \theta_0/\sqrt{1+\theta_0^2}}{rho = theta0/sqrt(1+theta0^2)}. The \code{get_cor} function returns the elements of the lower triangle of the correlation matrix, in column-major order.
##' @examples
##' th0 <- 0.5
##' stopifnot(all.equal(get_cor(th0),th0/sqrt(1+th0^2)))
##' set.seed(101)
##' C <- get_cor(rnorm(21), return_val = "mat")
##' ## test: round-trip
##' stopifnot(all.equal(get_cor(put_cor(C), return_val = "mat"), C))
##' @export
get_cor <- function(theta, return_val = c("vec", "mat")) {
    return_val <- match.arg(return_val)
    n <-  get_matdim(length(theta))
    R <- diag(n)
    R[upper.tri(R)] <- theta
    R[] <- crossprod(R) # R <- t(R) %*% R
    scale <- 1 / sqrt(diag(R))
    R[] <- scale * R * rep(scale, each = n) # R <- cov2cor(R)
    if (return_val == "mat") return(R)
    return(R[lower.tri(R)])
}

##' @rdname get_cor
##' @param C a correlation matrix
##' @param input_val input a vector of correlation values from the lower triangle ("vec"), or the full correlation matrix ("mat")? 
##' @export
put_cor <- function(C, input_val = c("mat", "vec")) {
    input_val <- match.arg(input_val)
    if (input_val == "vec") {
        ## construct matrix
        M <- diag(get_matdim(length(C)))
        M[lower.tri(M)] <- C
        M[upper.tri(M)] <- t(M)[upper.tri(M)]
        C <- M
    }
    cc2 <- chol(C)
    scale <- diag(cc2)
    cc2 <- cc2 %*% diag(1/scale)
    cc2[upper.tri(cc2)]
}

hasRandom <- function(x) {
    pl <- getParList(x)
    return(length(unlist(pl[grep("^theta",names(pl))]))>0)
}

##' retrieve parameters by name or index
##' @param parm parameter specifier
##' @param object fitted glmmTMB object
##' @param full include simple dispersion parameter?
##' @param include_nonest include mapped parameter indices for non-estimated (mapped or rank-deficient/dropped) parameters?
##' @noRd
getParms <- function(parm=NULL, object, full=FALSE, include_nonest = FALSE) {
    vv <- vcov(object, full=TRUE, include_nonest = include_nonest)
    sds <- sqrt(diag(vv))
    pnames <- names(sds) <- rownames(vv)       ## parameter names (user-facing)
    ee <- object$obj$env

    ## don't use object$obj$env$random; we want to keep "beta" vals, which may be
    ## counted as "random" if using REML
    drop_rand <- function(x) x[!x %in% c("b", "bzi", "bdisp")]
    if (!include_nonest) {
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
            parm <- grep("betadisp", intnames, invert=TRUE)
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
            parm <- grep("^betadisp", intnames)
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


## simplified version of glmmTMB package checking
##' @param this_pkg downstream package being tested
##' @param dep_pkg upstream package on which \code{this_pkg} depends
##' @param dep_type "ABI" or "package"
##' @param built_version a \code{numeric_version} object indicating what version of \code{dep_pkg} was used to  build \code{this_pkg}
##' @param warn (logical) warn if condition not met?
##' @importFrom utils packageVersion
##' @noRd
check_dep_version <- function(this_pkg = "glmmTMB",  dep_pkg = "TMB", dep_type = "package",
                              built_version = .TMB.build.version,
                              warn = TRUE) {
    ## FIXME: replace by TMB.Version() when available ?
    cur_version <- packageVersion(dep_pkg)
    result_ok <- cur_version == built_version
    if (!result_ok) {
        warning(
            sprintf("%s version mismatch: \n", dep_type),
            sprintf("%s was built with %s %s version %s\n",
                    this_pkg, dep_pkg, dep_type, built_version),
            sprintf("Current %s %s version is %s\n",
                    dep_pkg, dep_type, cur_version),
            sprintf("Please re-install %s from source ", this_pkg),
            "or restore original ",
            sQuote(dep_pkg), " package",
            " (see '?reinstalling' for more information)"
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
##' to specify that computations should be done in parallel.) To further
##' trace OpenMP settings, use \code{options(glmmTMB_openmp_debug = TRUE)}.
##' @seealso \code{\link[TMB]{benchmark}}, \code{\link{glmmTMBControl}}
##' @return \code{TRUE} or \code{FALSE} depending on availability of OpenMP,
##' @aliases openmp
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
#' @param update_gauss_disp update \code{betadisp} from variance to SD parameterization?
#' @export
up2date <- function(oldfit, update_gauss_disp = FALSE) {
  openmp(1)  ## non-parallel/make sure NOT grabbing all the threads!
  obj <- oldfit$obj
  ee <- obj$env

  if (isNullPointer(oldfit$obj$env$ADFun$ptr)) {

      pars <- c(grep("last\\.par", names(ee), value = TRUE), "par",
                "parfull")

      ## using ee$parList() rather than ee$parameters should help
      ##  with mapped parameter s... ??
      params <- ee$parList()

      if (length(ee$map) > 0) {
          for (n in names(ee$map)) {
              ee$parameters[[n]] <- params[[n]]
          }
      }
      
      ## change name of thetaf to psi
      if ("thetaf" %in% names(params)) {
          ee$parameters$psi <- params$thetaf
          ee$parameters$thetaf <- NULL
          pars <- c(grep("last\\.par", names(ee), value = TRUE),
                    "par")
          for (p in pars) {
              if (!is.null(nm <- names(ee[[p]]))) {
                  names(ee[[p]])[nm == "thetaf"] <- "psi"
              }
          }
      }
      if ("betad" %in% names(ee$parameters)) { #FIXME: DRY
      	ee$parameters$betadisp <- params$betad
      	ee$parameters$betad <- NULL
      	pars <- c(grep("last\\.par", names(ee), value = TRUE),
      						"par")
      	for (p in pars) {
      		if (!is.null(nm <- names(ee[[p]]))) {
      			names(ee[[p]])[nm == "betad"] <- "betadisp"
      		}
      	}
      	ee$data$Xdisp <- ee$data$Xd
      	ee$data$Xd <- NULL
      	ee$data$dispoffset <- ee$data$doffset
      	ee$data$doffset <- NULL
      }
      if(!"Zdisp" %in% names(ee$data)) {
      	ee$data$Zdisp <- new("dgTMatrix",Dim=c(as.integer(nrow(ee$data$Xdisp)),0L)) ## matrix(0, ncol=0, nrow=nobs)
      	ee$parameters$bdisp <- rep(0, ncol(ee$data$Zdisp))
      	ee$parameters$thetadisp <- numeric(0)
      }
      ee2 <- oldfit$sdr$env
      if ("thetaf" %in% names(ee2$parameters)) {
          ee2$parameters$psi <- ee2$parameters$thetaf
          ee2$parameters$thetaf <- NULL
      }

      for (i in seq_along(ee$data$terms)) {
          ee$data$terms[[i]]$simCode <- .valid_simcode[["random"]]
      }
      for (i in seq_along(ee$data$termszi)) {
          ee$data$termszi[[i]]$simCode <- .valid_simcode[["random"]]
      }
      
      if ("betad" %in% names(ee2$parameters)) { #FIXME: DRY
      	ee2$parameters$betadisp <- ee2$parameters$betad
      	ee2$parameters$betad <- NULL
      }

      ## prior_ivars, prior_fvars are defined in priors.R
      if (!"prior_distrib" %in% names(ee$data)) {
          ## these are DATA_IVECTOR but apparently after processing
          ##  TMB turns these into numeric ... ??
          for (v in prior_ivars) ee$data[[v]] <- numeric(0)
          for (v in prior_fvars) ee$data[[v]] <- numeric(0)

      }

      ## switch from variance to SD parameterization
      if (update_gauss_disp &&
          family(oldfit)$family == "gaussian") {
          ee$parameters$betadisp <- params$betadisp/2
          for (p in pars) {
              if (!is.null(nm <- names(ee[[p]]))) {
                  ee[[p]][nm == "betadisp"] <- ee[[p]][nm == "betadisp"]/2
              }
              if (!is.null(nm <- names(oldfit$fit[[p]]))) {
                  oldfit$fit[[p]][nm == "betadisp"] <- oldfit$fit[[p]][nm == "betadisp"]/2
              }
          }
      }

      oldfit$obj <- with(ee,
                       TMB::MakeADFun(data,
                                      parameters,
                                      map = map,
                                      random = random,
                                      silent = silent,
                                      DLL = "glmmTMB"))
      oldfit$obj$env$last.par.best <- ee$last.par.best
      ##
  }

  for (t in c("condReStruc", "ziRestruc", "dispRestruc")) {
      for (i in seq_along(oldfit$modelInfo$reStruc[[t]])) {
          oldfit$modelInfo$reStruc[[t]][[i]]$simCode <- .valid_simcode[["random"]]
      }
  }

  ## changed format of 'parallel' control to add autopar info
  if (length(p <- oldfit$modelInfo$parallel) <= 1) {
      if (!(is.null(p) || is.numeric(p))) {
          stop("oldfit$modelInfo$parallel has an unexpected value")
      }
      oldfit$modelInfo$parallel <- list(n = p, autopar = NULL)
  }

  ## dispersion was NULL rather than 1 in old R versions ...
  omf <- oldfit$modelInfo$family
  if (getRversion() >= "4.3.0" &&
      !("dispersion" %in% names(omf))) {
      ## don't append() or c(), don't want to lose class info
      oldfit$modelInfo$family$dispersion <- 1
  }
  if (!"priors" %in% names(oldfit$modelInfo)) {
      ## https://stackoverflow.com/questions/7944809/assigning-null-to-a-list-element-in-r
      ## n.b. can't use ...$priors <- NULL
      oldfit$modelInfo["priors"] <- list(NULL)
  }

  if ("Xd" %in% names(ee$data)) {
      ee$data[["Xdisp"]] <- ee$data[["Xd"]]
      ee$data[["Xd"]] <- NULL
  }

  if ("XdS" %in% names(ee$data)) {
      ee$data[["XdispS"]] <- ee$data[["XdS"]]
      ee$data[["XdS"]] <- NULL
  }

  return(oldfit)
}

#' Load data from system file, updating glmmTMB objects
#'
#' @param fn partial path to system file (e.g. test_data/foo.rda)
#' @param verbose print names of updated objects?
#' @param mustWork fail if file not found?
#' @param \dots values passed through to \code{up2date}
#' @export
gt_load <- function(fn, verbose=FALSE, mustWork = FALSE, ...) {
    sf <- system.file(fn, package = "glmmTMB")
    found_file <- file.exists(sf)
    if (mustWork && !found_file) {
        stop("couldn't find system file ", sf)
    }

    L <- load(sf)
    for (m in L) {
        if (inherits(get(m), "glmmTMB")) {
            if (verbose) cat(m,"\n")
            assign(m, up2date(get(m), ...))
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
par_components <- c("beta","betazi","betadisp","theta","thetazi","psi")

## all parameters, including both mapped and rank-dropped
getParnames <- function(object, full, include_dropped = TRUE, include_mapped = TRUE) {
                           
  mkNames <- function(tag="") {
      X <- getME(object,paste0("X",tag))
      if (is.null(X) && tag == "disp") stop("are you using a stored model from an earlier glmmTMB version? try running up2date()")
      dropped <- attr(X, "col.dropped") %||% numeric(0)
      ntot <- ncol(X) + length(dropped)
      ## identical instead of ==; ncol(X) may be NULL for older models
      if (ntot == ncol(X) || !include_dropped) {
          nn <- colnames(X)
      } else {
          nn <- character(ntot)
          nn[-dropped] <- colnames(X)
          nn[ dropped] <- names(dropped)
      }
      if (trivialFixef(nn, tag)
          ## if 'full', keep disp even if trivial, if used by family
          && !(full && tag =="disp" &&
               (usesDispersion(family(object)$family) && !zeroDisp(object)))) {
          return(character(0))
      }
      if (tag == "") return(nn)
      return(paste(tag,nn,sep="~"))
  }

  nameList <- setNames(Map(mkNames, c("", "zi", "disp")),
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
                    list(theta = reNames("cond"), 
                    		 thetazi = reNames("zi"), 
                    		 thetadisp = reNames("disp")))

      ##
      if (length(fp <- family_params(object)) > 0) {
          nameList <- c(nameList, list(psi = names(fp)))
      }
      
  }

  if (!include_mapped) {
     map <- object$obj$env$map
     if (length(map)>0) {
         for (m in seq_along(map)) {
            if (length(NAmap <- which(is.na(map[[m]])))>0) {
                w <- match(names(map)[m], par_components) ##
                if (length(nameList)>=w) { ## may not exist if !full
                    nameList[[w]] <- nameList[[w]][-NAmap]
                }
            }
         } ## for (m in seq_along(map))
     } ## if (length(map) > 0)
  }

  return(nameList)
}

## OBSOLETE (?)

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

collapse_list <- function(pList) {
    ## workaround to get non-unique names ...
    pList <- mapply(function(x, n) { setNames(x, rep(n, length(x))) },
                    pList, names(pList))
    pvec <- unlist(unname(pList))
}

make_pars <- function(pList, ..., include_extra = TRUE) {
    ## FIXME: check for name matches, length matches etc.
    ## (useful errors)
    ## better to split by name first??

    L <- list(...)
    unmatched <- setdiff(names(L), names(pList))
    if (length(unmatched) > 0) {
        warning(sprintf("unmatched parameter names: %s",
                        paste(unmatched, collapse =", ")))
    }

    if (!include_extra) L <- L[intersect(names(L), names(pList))]
    for (nm in names(L)) {
        if ((len1 <- length(pList[[nm]])) == (len2 <- length(L[[nm]]))) {
            pList[[nm]] <- L[[nm]]
        } else {
            ## skip cases with different length (== partially-mapped vectors)
            plen <- len1 + length(attr(pList[[nm]], "map"))
            ## if L[[nm]] is a list, that's because we're
            ##  passing a partial b vector ...
            if (plen != len2 && !is.list(L[[nm]])) {
                warning(sprintf("length mismatch in component %s (%d != %d); not setting",
                                nm, len1, len2))
            }
        }
    }
    return(collapse_list(pList))
}

##' helper function to modify simulation settings for random effects
##'
##' This modifies the TMB object \emph{in place} (beware!)
##' Ultimately this will allow \code{terms} to be a vector of term names,
##' with a matching \code{val} vector to specify the behaviour for each term
##'
##' @param g a TMB object
##' @param val a legal setting for sim codes ("zero", "random", or "fix")
##' @param terms which terms to apply this to
##' @export
##' 
set_simcodes <- function(g, val = "zero", terms = "ALL") {
    ee <- g$env
    if (terms != "ALL") stop("termwise setting of simcodes not implemented yet")
    if (terms == "ALL") {
        for (i in seq_along(ee$data$terms)) {
            ee$data$terms[[i]]$simCode <- .valid_simcode[[val]]
        }
    }

}

##' Simulate from covariate/metadata in the absence of a real data set (EXPERIMENTAL)
##'
##' See \code{vignette("sim", package = "glmmTMB")} for more details and examples,
##' and \code{vignette("covstruct", package = "glmmTMB")}
##' for more information on the parameterization of different covariance structures.
##'
##' @inheritParams glmmTMB
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
##' (\code{beta}, \code{betazi}, \code{betadisp}, \code{theta}, etc.) to
##' be used in the model. If \code{b} is specified in this list, then the conditional modes/BLUPs
##' will be set to these values; otherwise they will be drawn from the appropriate Normal distribution
##' @param ... other arguments to \code{glmmTMB} (e.g. \code{family})
##' @param return_val what information to return: "sim" (the default) returns a list of vectors of simulated outcomes; "pars" returns the default parameter vector (this variant does not require \code{newparams} to be specified, and is useful for figuring out the appropriate dimensions of the different parameter vectors); "object" returns a fake \code{glmmTMB} object (useful, e.g., for retrieving the Z matrix (\code{getME(simulate_new(...), "Z")}) or covariance matrices (\code{VarCorr(simulate_new(...))}) implied by a particular set of input data and parameter values)
##' @examples
##' ## use Salamanders data for structure/covariates
##' sim_count <- simulate_new(~ mined + (1|site),
##'              newdata = Salamanders,
##'              zi = ~ mined,
##'              family = nbinom2,
##'              newparams = list(beta = c(2, 1),
##'                          betazi = c(-0.5, 0.5), ## logit-linear model for zi
##'                          betadisp = log(2), ## log(NB dispersion)
##'                          theta = log(1)) ## log(among-site SD)
##' )
##' sim_obj <- simulate_new(~ mined + (1|site),
##'             return_val = "object",
##'              newdata = Salamanders,
##'              zi = ~ mined,
##'              family = nbinom2,
##'              newparams = list(beta = c(2, 1),
##'                          betazi = c(-0.5, 0.5), ## logit-linear model for zi
##'                          betad = log(2), ## log(NB dispersion)
##'                          theta = log(1)) ## log(among-site SD)
##' )
##' data("sleepstudy", package = "lme4")
##' sim_obj <- simulate_new(~ 1 + (1|Subject) + ar1(0 + factor(Days)|Subject),
##'              return_val = "pars",
##'              newdata = sleepstudy,
##'              family = gaussian,
##'              newparams = list(beta = c(280, 1),
##'                          betad = log(2), ## log(SD)
##'                          theta = log(c(2, 2, 1))),
##' )
##' 
##' @export
simulate_new <- function(object,
                         nsim = 1,
                         seed = NULL,
                         family = gaussian,
                         newdata, newparams, ...,
                         return_val = c("sim", "pars", "object")) {
    return_val <- match.arg(return_val)
    family <- get_family(family)
    ## truncate
    if (length(object) == 3) stop("simulate_new should take a one-sided formula")
    newparams0 <- newparams
    ## store original params
    ## (in case we need both complete-b and unmapped-b versions)
    
    ## fill in fake LHS
    form <- object
    form[[3]] <- form[[2]]
    form[[2]] <- quote(..y)
    ## insert a legal value: 1.0 is OK as long as family != "beta"
    ## (note the family *function* is 'beta_family' but the internal
    ##  $family value is 'beta')
    newdata[["..y"]] <- if (family$family == "beta") 0.5 else 1.0
    r1 <- glmmTMB(form,
                  data = newdata,
                  family = family,
                  ## make sure optim doesn't actually do anything
                  ## (if return_val is "object")
                  control = glmmTMBControl(optCtrl = list(iter.max = 0)),
                  ...,
                  doFit = FALSE)
    ## sort out components of b (if necessary)
    if ("b" %in% names(newparams)) {
        components <- c("cond", "zi")
        cnames <- paste0(components, "ReStruc")
        tnames <- c("terms", "termszi")
        if (!is.list(newparams$b)) {
            b_inds <- seq_along(newparams$b)
            for (i in seq_along(cnames)) { ## loop over RE categories in object ('cond' and 'zi')
                for (j in seq_along(r1[[cnames[i]]])) { ## loop over terms within the category
                    r1[[cnames[i]]][[j]]$simCode <- .valid_simcode[["fix"]]
                    r1$data.tmb[[tnames[i]]][[j]]$simCode <- .valid_simcode[["fix"]]
                }
            }
        } else {
            restrucs <- r1[cnames]
            b_inds <- get_b_inds(restrucs, names(newparams$b))
            b_terms <- get_b_inds(restrucs, names(newparams$b),
                                  ret_val = "terms")
            ##  b_terms gives indices in overall sequence, not
            ## indices per term -- we have to split these by
            ## term. Should probably be handled upstream?
            for (i in seq_along(cnames)) {
                nre <- length(r1[[i]])
                cur_b <- which(b_terms < nre)
                if (length(cur_b) > 0) {
                    for (j in b_terms[cur_b]) {
                        ## ugh, have to set in two places: do this upstream
                        ## of glmmTMB() call?
                        r1[[cnames[i]]][[j]]$simCode <- .valid_simcode[["fix"]]
                        r1$data.tmb[[tnames[i]]][[j]]$simCode <- .valid_simcode[["fix"]]
                    }
                    b_terms <- b_terms[-cur_b]
                }
                b_terms <- b_terms - nre
            }
        }
        n_b <- length(r1$parameters$b)
        b_fac <- seq(n_b)
        b_fac[unlist(b_inds)] <- NA
        ## TMB complains if number of levels doesn't match values:
        ##  make into a factor *after* setting NA values
        b_fac <- factor(b_fac)
        new_b <- rep(0, n_b)
        new_b[unlist(b_inds)] <- unlist(newparams$b)
        r1$parameters$b <- new_b
        r1$map <- r1$mapArg <- list(b = b_fac)
    }
    ## construct TMB object, but don't fit it
    ## (for cnms etc., simulations, etc.)
    r2 <- fitTMB(r1, doOptim = FALSE)

    set_b <-  function(x, b, map = NULL) {
        if (!is.null(map)) {
            b <- b[!is.na(map)]
        }
        pList <- split(x, names(x))
        pList$b <- b
        return(collapse_list(pList))
    }

    pars <- do.call("make_pars",
                    c(list(r2$env$parameters), newparams,
                      list(include_extra = FALSE)))

    if (!is.null(seed)) set.seed(seed)
    if (return_val %in% c("pars", "object")) {
        b_vals <- r2$simulate(par = pars)$b
        if (return_val == "pars") {
            pars <- do.call("make_pars",
                            c(list(r2$env$parameters), newparams,
                              list(include_extra = TRUE)))
            return(set_b(pars, b_vals))
        }
        ## FIXME: need to set start= values as well when simulating
        ## object??
        ## what do we do when b is specified/partially specified?
        for (nm in names(newparams)) {
            r1$parameters[[nm]] <- newparams[[nm]]
        }
        r3 <- suppressWarnings(fitTMB(r1, doOptim = TRUE))
        r3$fit$parfull <- set_b(r3$fit$parfull, b_vals, r1$map$b)
        return(r3)
    }
    replicate(nsim, r2$simulate(par = pars)$yobs, simplify = FALSE)
}

set_class <- function(x, cls, prepend = TRUE) {
    if (is.null(x)) return(NULL)
    if (!prepend) class(x) <- cls
    else class(x) <- c(cls, class(x))
    x
}

## convert from parameter name to component name or vice versa
## first name shoudl be em
compsyn <- c(cond = "", zi = "zi", disp = "d")
match_names <- function(x, to_parvec = FALSE, prefix = "beta") {
    if (to_parvec) {
        ## "cond" -> "theta" etc.
        return(paste0(prefix, compsyn[x]))
    } else {
        ## "beta" -> "cond" etc.
        x <- gsub(prefix, "", x)
        return(names(compsyn)[match(x, compsyn)])
    }
}

get_family <- function(family) {
    if (is.character(family)) {
        if (family=="beta") {
            family <- "beta_family"
            warning("please use ",sQuote("beta_family()")," rather than ",
                    sQuote("\"beta\"")," to specify a Beta-distributed response")
        }
        family <- get(family, mode = "function", envir = parent.frame(2))
    }

    if (is.function(family)) {
        ## call family with no arguments
        family <- family()
    }

    ## FIXME: what is this doing? call to a function that's not really
    ##  a family creation function?
    if (is.null(family$family)) {
      print(family)
      stop("after evaluation, 'family' must have a '$family' element")
    }
    return(family)
}

## ugh, better way to do this?
get_re_names <- function(re) {
    names(unlist(sapply(re, function(x) sapply(x, function(y) NA))))
}

#' @param reStrucs a list containing conditional and z-i RE structures
#' @param b_names vector of names matching RE terms
#' @examples
#' data("sleepstudy", package = "lme4")
#' fm1 <- glmmTMB(Reaction ~ 1 + (1|Subject) + ar1(0+factor(Days)|Subject), sleepstudy)
#' re <- fm1$modelInfo$reStruc
#' get_b_inds(re, "1|Subject")
#' @noRd
get_b_terms <- function(nms, inms) {
    squash_ws <- function(x) gsub(" ", "", x)
    nms <- squash_ws(nms)
    inms <- squash_ws(inms)
    w <- match(nms, inms)
    unmatched <- which(is.na(w))
    if (length(unmatched)>0) {
        w[unmatched] <- match(nms[unmatched],
                              gsub("^(cond|zi)ReStruc\\.", "", inms))
    }
    if (any(is.na(w))) {
        stop("unmatched RE terms")
    }
    return(w)
}
                        
get_b_inds <- function(reStrucs, b_names, ret_val = c("indices", "terms")) {
    ret_val <- match.arg(ret_val)
    bfun <- function(x) {
        if (length(x) == 0) return(numeric(0))
        with(x, blockSize * blockReps)
    }
    retrms <- unlist(sapply(reStrucs, function(x) sapply(x, bfun)))
    ## set up indices ...
    if (ret_val == "terms") {
        return(get_b_terms(b_names, get_re_names(reStrucs)))
    }
    w <- get_b_terms(b_names, get_re_names(reStrucs))
    inds <- cumsum(c("start" = 0, retrms))
    ## first try to match full name (component + term)
    ## set specified values
    res <- lapply(w, function(i) seq(inds[i]+1, inds[i+1]))
    names(res) <- b_names
    res
}

## add negative-value check to binomial initialization method
our_binom_initialize <- function(family) {
    newtest <- substitute(
        ## added test for glmmTMB
        if (any(y<0)) {
            stop(sprintf('negative values not allowed for the %s family', FAMILY))
        }
      , list(FAMILY=family))
    b0 <- binomial()$initialize
    b0[[length(b0)+1]] <- newtest
    return(b0)
}

