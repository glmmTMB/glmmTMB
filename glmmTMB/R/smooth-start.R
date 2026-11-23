## Scale-aware starts for conditional smooths in the smooth2random basis.
## Keep this separate from startParams(): explicit starts are applied afterwards.
smoothStart <- function(parameters, condList, condReStruc, data.tmb,
                        family, map = NULL) {
    info <- condList$reTrms$smooth_info
    positions <- which(lengths(info) > 0L)
    if (!length(positions) || !is.null(map$theta)) return(parameters)
    links <- switch(family$family,
                    gaussian = "identity",
                    poisson =, nbinom1 =, nbinom2 =, tweedie = "log",
                    binomial =, beta = "logit",
                    Gamma = c("log", "inverse"), character())
    if (!(family$link %in% links)) return(parameters)

    ## Retain the original scope for Gaussian models. Non-Gaussian models
    ## can have zero-inflation or dispersion predictors, whose starts are
    ## left alone: only the conditional smooth parameters are initialized.
    if (family$family == "gaussian") {
        if (ncol(data.tmb$Xzi) > 0L || ncol(data.tmb$XziS) > 0L ||
            ncol(data.tmb$Zzi) > 0L || length(parameters$betadisp) != 1L ||
            ncol(data.tmb$Zdisp) > 0L || !is.null(map$betadisp)) return(parameters)
        Xdisp <- if (ncol(data.tmb$Xdisp)) data.tmb$Xdisp else data.tmb$XdispS
        if (any(Xdisp != 1) || any(data.tmb$dispoffset != 0)) return(parameters)
    }

    smooths <- lapply(info[positions], `[[`, "sm")
    ## Each supported smooth block has one penalty and one log SD. Prediction
    ## reconstruction has no sm object and must not initialize from new data.
    if (any(vapply(smooths, function(sm)
        is.null(sm$X) || length(sm$S) != 1L, logical(1))) ||
        any(vapply(condReStruc[positions], function(block)
            block$blockNumTheta != 1L, logical(1)))) return(parameters)

    y <- data.tmb$yobs
    ## Only identity-link offsets can be subtracted on the response scale.
    if (family$family == "gaussian") y <- y - data.tmb$offset
    w <- data.tmb$weights
    keep <- is.finite(y) & is.finite(w) & w > 0
    if (sum(keep) < 2L) return(parameters)
    y <- y[keep]
    w <- w[keep]
    ## This is a scale heuristic, not a family-specific dispersion estimate.
    ## It only scales smooth theta: the dispersion start stays unchanged.
    response_var <- sum(w * (y - stats::weighted.mean(y, w))^2) / sum(w)
    if (!is.finite(response_var) || response_var <= 0) return(parameters)

    ## Use the original constrained bases and their matching penalties.
    ## The penalty is already scaled by smoothCon(); smooth2random() makes
    ## its penalized part an identity matrix. For Gaussian models this is
    ## motivated by sigma_b^2 = phi / lambda; for GLMs it is an initial guess.
    X <- do.call(cbind, lapply(smooths, function(sm) sm$X[keep, , drop = FALSE]))
    S <- lapply(smooths, function(sm) sm$S[[1L]])
    off <- head(cumsum(c(1L, vapply(S, ncol, integer(1)))), -1L)
    lambda <- tryCatch(mgcv::initial.sp(sqrt(w) * X, S, off),
                       error = function(e) NULL)
    if (is.null(lambda)) return(parameters)
    theta <- (log(response_var) - log(lambda)) / 2
    if (any(!is.finite(theta))) return(parameters)

    ## Random-effect blocks can have several covariance parameters and need
    ## not precede smooths in the formula. Only change the smooth positions.
    sizes <- vapply(condReStruc, `[[`, numeric(1), "blockNumTheta")
    theta_pos <- cumsum(c(0L, sizes))[positions] + 1L
    parameters$theta[theta_pos] <- theta
    parameters
}
