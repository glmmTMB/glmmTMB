#' @param model fitted TMB model
#' @param data original data
#' @param epsilon perturbation size
#' @param inds indices of observations to perturb
#' @param fit_method "hack" = quick, without finalizing model object; "update" = slow, complete model fit + update
#' @param pred_method "hack" = based on linear model; "predict" = from fitted model object
#' @param scale data ("response") or linear predictor ("link") scale
#' @param progress progress bar?
#' @param opt_args additional arguments for optimizer
#' @param return_delta for diagnostics/debugging: return unscaled delta rather than delta/eps?
leverage_brute <- function(model, data, epsilon = 1e-3, inds = seq(nrow(data)),
                           fit_method = c("hack", "update"),
                           pred_method = c("hack", "predict"),
                           scale = c("response", "link"),
                           progress = FALSE,
                           opt_args = list(),
                           return_delta = FALSE
                           ) {

    scale <- match.arg(scale)
    fit_method <- match.arg(fit_method)
    pred_method <- match.arg(pred_method)

    n <- length(inds)

    if (progress) pb <- txtProgressBar(max = n, style = 3)
    ## for now, compute all leverages on link scale
    y_pred <- predict(model, type = "link")
    leverage <- rep(NA_real_, n)
    yname  <- as.character(formula(model)[[2]])

    ## extract parameters so we can start from best vals
    p0 <- with(model$obj$env, parList(last.par.best[-random]))
    p0 <- p0[lengths(p0) > 0]
    p0 <- p0[setdiff(names(p0), "b")]  ## drop 'b' parameters

    X <- getME(model, "X")
    Z <- getME(model, "Z")

    for (j in seq_along(inds)) {
        if (progress) setTxtProgressBar(pb, j)
        i <- inds[j]
        data_perturb <- data
        data_perturb[[yname]][i] <-  data[[yname]][i] + epsilon
        if (fit_method == "hack") {
            ## quick/hacked new fit
            ## don't want to see warnings about non-integer counts
            suppressWarnings(
                newfit0 <- update(model, start = p0, data = data_perturb, verbose = FALSE, doFit = FALSE)
            )
            system.time(newfit1 <- fitTMB(newfit0, doOptim = FALSE)) ## 1 second
            system.time(newfit2 <- with(newfit1,
                                        do.call(nlminb,
                                                c(list(start=par, objective=fn, gradient=gr),
                                                  opt_args)))
                        )
        } else {
            ## full new fit
            suppressWarnings(
                newfit0 <- update(model, start = p0, data = data_perturb, verbose = FALSE)
            )
            newfit1 <- newfit0$obj
        }
        pp <- with(newfit1$env, parList(last.par.best[-random]))
        pp$b  <-  newfit1$report()$b
        if (pred_method == "hack") {
            y_pred_pert <- drop(X[i,] %*% pp[["beta"]] + Z[i,] %*% pp[["b"]])
        } else {
            if (fit_method == "hack") stop("can't do regular pred with hacked fit")
            y_pred_pert <- predict(newfit0, type = "link")[i]
        }
        if (scale == "response") {
            linkinv <- family(model)$linkinv
            y_pred_pert <- linkinv(y_pred_pert)
            y_pred[i] <- linkinv(y_pred[i])
        }
        leverage[j] <- (y_pred_pert - y_pred[i])
        if (!return_delta) leverage[j] <-  leverage[j] / epsilon
    }
    if (progress) close(pb)
    return(leverage)
}

## need to modify src/Makevars in glmmTMB directory to contain this
## (no fopenmp!)
install_glmmTMB <- function(pkgdir, libdir, clean_src = TRUE) {
    ## save existing src/Makevars, overwrite with what we want
    flags <- c("PKG_CPPFLAGS = -DTMBAD_FRAMEWORK -DTMBAD_INDEX_TYPE=uint64_t -DTMB_MAX_ORDER=4",
      "## PKG_LIBS = $(SHLIB_OPENMP_CXXFLAGS)",
      "## PKG_CXXFLAGS=$(SHLIB_OPENMP_CXXFLAGS)")
    td <- tempdir()
    makevars <- file.path(pkgdir, "src", "Makevars")
    file.rename(makevars, file.path(td, "Makevars"))
    on.exit(file.rename(file.path(td, "Makevars"), makevars))
    unlink(makevars)
    writeLines(flags, makevars)
    if (clean_src) {
        unlink(list.files(file.path(pkgdir, "src"),
                          pattern="\\.(o|so)$",
                          full.names = TRUE))
    }
    if (!dir.exists(libdir)) dir.create(libdir)
    system(sprintf("R CMD INSTALL -l %s %s",
                   libdir, pkgdir))
}

## need to install this way (with location adjusted to your liking)
## R CMD INSTALL -l ~/students/agronah/reduced_rank_mm/glmmTMB_lev glmmTMB

## from https://github.com/glmmTMB/glmmTMB/blob/leverage/misc/leverage.R
## @param diag Get diagonal only?
leverage <- function(fm, diag=TRUE) {
    library(RTMB)
    library(Matrix)
    has.random <- any(fm$obj$env$lrandom())
    obj <- fm$obj
    ## We mess with these... (cleanup on exit!)
    restore.on.exit <- c("ADreport",
                         "parameters",
                         "data")
    oldvars <- sapply(restore.on.exit, get, envir=obj$env, simplify=FALSE)
    restore.oldvars <- function(){
        for(var in names(oldvars)) assign(var, oldvars[[var]], envir=obj$env)
    }
    on.exit({restore.oldvars(); obj$retape()})
    ## #################################################################
    ## Get derivatives of prediction
    ##
    ##    mu_hat( b_hat( theta_hat(yobs), yobs) , theta_hat(yobs) )
    ##
    ## wrt. yobs.
    ## Note the three partial derivative 'paths' to consider.
    ## We can get this derivative by
    ##  1. yobs -> theta_hat(yobs)
    ##  2. theta -> mu_hat( b_hat( theta, yobs) , theta ) [fixed yobs]
    ##  3. yobs -> mu_hat( b_hat( theta, yobs) , theta ) [fixed theta]
    ## #################################################################
    ##parhat <- obj$env$last.par.best
    parhat <- fm$fit$parfull
    pl <- obj$env$parList(par=parhat)
    yobs <- obj$env$data$yobs
    obj$env$parameters <- pl
    theta <- parhat[!obj$env$lrandom()]  ## ALL top-level parameters
    b <- parhat[obj$env$lrandom()]
    if (!is.null(obj$env$spHess)) {
        Hbb <- obj$env$spHess(parhat, random=TRUE) ## Needed later for RE models
    }
    ## #################################################################
    ## 1. Get partial derivatives of theta_hat wrt to yobs
    ## Note: length(yobs) much greater that length(theta)
    ##       ==> Reverse mode AD is suitable !
    ## #################################################################
    ## Move 'yobs' from data -> parameters (preserve C++ template order!)
    obj$env$parameters <- c(list(yobs = yobs), obj$env$parameters)
    obj$env$data$yobs <- NULL
    obj$retape()
    ## New objective parameters: (yobs, b, theta)
    nobs <- length(obj$env$parameters$yobs)
    nb <- length(obj$env$random)
    ntheta <- length(obj$env$par) - nobs - nb
    TMB::config(tmbad.atomic_sparse_log_determinant=0, DLL="RTMB") ## TMB FIXME
    F <- GetTape(obj)
    r <- obj$env$random ## Which are random
    p <- tail(1:(nobs+ntheta), ntheta) ## Which are parameters *after* removing random
    ThetaHat <- F$laplace(r)$newton(p)
    J <- ThetaHat$jacobian(ThetaHat$par())
    ## Extra stuff we need in (3)
    F. <- F$jacfun() ## (yobs, [b], theta) -> (yobs, [b], theta)
    F. <- MakeTape(function(y) F.( c(y, parhat) ) [r] , yobs) ## yobs -> b
    Hby <- F.$jacfun(sparse=TRUE)(yobs)
    ## #################################################################
    ## 2. Get partial derivatives of mu_hat wrt to theta for *fixed* yobs
    ## Note: length(mu) much greater that length(theta)
    ##       ==> Forward mode AD is suitable !
    ## #################################################################
    obj$env$data$yobs <- yobs
    obj$env$parameters$yobs <- NULL
    obj$retape()
    r <- obj$env$random ## Which are now random
    F <- GetTape(obj)
    Bhat <- F$newton(r) ## theta -> bhat
    obj$env$data$doPredict <- as.double(1) ## Enable prediction of 'mu'
    obj$env$data$whichPredict <- as.double(1:nobs)
    obj$env$ADreport <- TRUE ## Modify return value from Cpp
    obj$retape()
    F <- GetTape(obj) ## (b, theta) -> mu
    ## This doesn't work:
    ## MuHat <- MakeTape(function(theta)F(c(Bhat(theta), theta)), theta)
    MuHat <- MakeTape(function(theta) {
        par <- advector(nb + ntheta) ## glmmTMB mixes order of parameters and random effects...
        r <- obj$env$lrandom()
        par[r] <- Bhat(theta)
        par[!r] <- theta
        F(par)
    } , theta)
    ## 'Adjoint trick'
    T2 <- MakeTape(function(weight) {
        WT <- MakeTape(function(th) sum(MuHat(th) * weight), theta)
        WT$jacfun()(advector(theta))
    }, rep(1,nobs))
    J2 <- T2$jacobian(yobs)
    if (diag) {
        term1 <- colSums(J*J2)
    } else {
        term1 <- t(J) %*% J2
    }
    ## #################################################################
    ## 3. Get partial derivatives of mu_hat wrt yobs for fixed theta
    ## Note: Tricky!
    ##       
    ## #################################################################
    term2 <- 0
    if (has.random) {
        F2 <- MakeTape(function(b) {
            par <- advector(nb + ntheta) ## glmmTMB mixes order of parameters and random effects...
            r <- obj$env$lrandom()
            par[r] <- b
            par[!r] <- theta
            F(par)
        }, b) ## (b) -> mu
        Hmb <- F2$jacfun(sparse=TRUE)(b) ## sparse deriv mu wrt b
        ## Implicit function theorem gives final partial deriv matrix:
        ##   - Hby %*% solve(Hbb) %*% Hbm
        ## of which we need the diagonal.
        ## Because mu and yobs link to the same random effects, all required b-cliques are part of Hbb !
        ## It follows that we can replace solve(Hbb) by its subset iH !
        if (diag) {
            if (length(b) == 0) {
                iH <- solve(Hbb)
            } else {
                iH <- TMB:::solveSubset(Hbb)
            }
            term2 <- -colSums( Hby * (  iH %*% t(Hmb) ) )
        } else {
            term2 <- -t(Hby) %*% solve(Hbb, t(Hmb))
        }
    }
    term1 + term2
}

peakRAM_testfun <- function(nsubj = 100, ntax = 100, d = 2,
                            include_ttt = FALSE, seed = 101,
                            vars_include = c("nsubj", "ntax", "d", "include_ttt")) {

    if (packageVersion("glmmTMB") < "1.1.11") stop("are you using the formula_env branch?")
    ## must assume we are using the `formula_env` branch of glmmTMB!! hard to test though

    library(peakRAM)

    ##
    set.seed(seed)
    dd <- expand.grid(subject = factor(seq(nsubj)),
                      taxon = factor(seq(ntax)))
    dd$group <- factor(ifelse(as.numeric(dd$subject) < nsubj %/% 2, "a", "b"))

    if (include_ttt) {
        form <- y  ~ 1 + us(1 + group | taxon) + rr(0 + taxon | subject, d)
    } else {
        form <- y ~ 1 + rr(0 + taxon | subject, d)
    }
    dd$y <- simulate_new( form[-2],
                 family = nbinom2,
                 newdata = dd,
                 control = list(set_formula_env = FALSE),
                 newparams = list(beta = 1,
                                  betadisp = 1,
                                  theta = rep(0.1, ntheta(nsubj, ntax, d, include_ttt)))
                 )[[1]]
    ## have to run this without parallelization, since we had to turn off
    ## OpenMP for leverage calculations (we could try to load the full
    ## version of glmmTMB with autopar for fitting the model, then
    ## detach and load the glmmTMB_lev
    p1 <- peakRAM(
        mod <- glmmTMB(form,
                       family = nbinom2,
                       data = dd)
    )
    tmpf <- function(x, task = "model_fit") {
        names(x) <- c("task", "time_sec", "total_RAM_Mb", "peak_RAM_Mb")
        x$task <- task
        v <- mget(vars_include, inherits = TRUE)
        x <- do.call(data.frame, c(v, list(x)))
        return(x)
    }
    p2 <- peakRAM(leverage(mod))
    rbind(tmpf(p1), tmpf(p2, task = "leverage"))
}


ntheta <- function(nsubj = 100, ntax = 100, d = 2, include_ttt = FALSE) {
    rr_n <- ntax*d - choose(d,2)
    rr_n + ifelse(include_ttt, 3, 0)
}
