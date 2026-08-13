#' compute denominator degrees-of-freedom approximations
#'
#' \code{dof_KR} uses an adaptation of the machinery from the \code{pbkrtest} package
#' to compute the Kenward-Roger approximation of the 'denominator degrees of freedom' for
#' each fixed-effect coefficient in the conditional model; \code{dof_satt} does the same
#' for Satterthwaite approximations
#' @return a named vector of ddf for each conditional fixed-effect parameter; \code{dof_KR} includes attributes 'vcov'
#' (Kenward-Roger adjusted covariance matrix) and 'se' (the corresponding standard errors)
#' @details Kenward-Roger adjustments \emph{should not be used} for models fitted with ML rather than REML;
#' the theory is only well understood, and the model is only tested, for LMMs (\code{family = "gaussian"}).
#' Use at your own risk for GLMMs!
#' @param model a fitted \code{glmmTMB} object
#' @export
## avoid conflict with insight::dof_kenward ...
## FIXME: check with various combinations of mapping etc.
dof_KR <- function(model) {
    fe <- fixef(model)$cond
    param_names <- names(fe)
    L <- as.data.frame(diag(rep(1, length(fe))))
    krvcov <- .vcov_kenward_adjusted(model)

    dof <- vapply(L, .kenward_adjusted_ddf, model = model, adjusted_vcov = krvcov,
                  FUN.VALUE = numeric(1))
    names(dof) <- param_names
    attr(dof, "vcov") <- krvcov
    attr(dof, "se") <- abs(sqrt(diag(krvcov)))
    dof
}

## The following code was taken from the "pbkrtest" package and slightly modified
## and then slightly modified again for "glmmTMB"
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}

## FIXME: use built-in family$var functions? Does sigma^2*family()$variance() work universally?
.family_var_func<-function(model){
    if(model$modelInfo$family$family=="gaussian"){
        (predict(model, type = "disp")^2)
    }else if(model$modelInfo$family$family=="nbinom1"){
        (predict(model, type="response")*(1+predict(model, type="disp")))
    }else if(model$modelInfo$family$family=="nbinom2"){
        (predict(model, type="conditional")+(predict(model, type="conditional")^2/predict(model, type="disp")))
    }else if(model$modelInfo$family$family=="nbinom12"){
        (predict(model, type="conditional")*(1+predict(model, type="disp")+(predict(model, type="conditional")/.link_to_response(model$fit$par["psi"],model))))
    }else if(model$modelInfo$family$family=="Gamma"){
        (predict(model, type="response")*predict(model, type="disp"))
    }
}

## FIXME: use family()$linkfun
.link_func<-function(x, model){
    if(model$modelInfo$family$link=="identity"){
        x
    }else if(model$modelInfo$family$link=="log"){
        log(x)
    }else if(model$modelInfo$family$link=="inverse"){
        1/x
    }
}

## FIXME: use family()$linkinv
.link_to_response<-function(x,model){
    if(model$modelInfo$family$link=="identity"){
        x
    }else if(model$modelInfo$family$link=="log"){
        exp(x)
    }else if(model$modelInfo$family$link=="inverse"){
        1/x
    }
}

.kenward_adjusted_ddf <- function(model, linear_coef, adjusted_vcov) {
    .adjusted_ddf(adjusted_vcov, linear_coef, stats::vcov(model)$cond)
}

.adjusted_ddf <- function(adjusted_vcov, linear_coef, unadjusted_vcov = adjusted_vcov) {

    if (!is.matrix(linear_coef)) {
        linear_coef <- matrix(linear_coef, ncol = 1)
    }
    vlb <- sum(linear_coef * (unadjusted_vcov %*% linear_coef))
    theta <- Matrix::Matrix(as.numeric(outer(linear_coef, linear_coef) / vlb), nrow = length(linear_coef))
    P <- attr(adjusted_vcov, "P")
    W <- attr(adjusted_vcov, "W")

    A1 <- A2 <- 0
    theta_unadjusted_vcov <- theta %*% unadjusted_vcov
    n.ggamma <- length(P)
    for (ii in 1:n.ggamma) {
        for (jj in ii:n.ggamma) {
            if (ii == jj) {
                e <- 1
            } else {
                e <- 2
            }
            ui <- as.matrix(theta_unadjusted_vcov %*% P[[ii]] %*% unadjusted_vcov)
            uj <- as.matrix(theta_unadjusted_vcov %*% P[[jj]] %*% unadjusted_vcov)
            A1 <- A1 + e * W[ii, jj] * (sum(diag(ui)) * sum(diag(uj)))
            A2 <- A2 + e * W[ii, jj] * sum(ui * t(uj))
        }
    }

    B <- (A1 + 6 * A2) / 2
    g <- (2 * A1 - 5 * A2) / (3 * A2)
    c1 <- g / (3 + 2 * (1 - g))
    c2 <- (1 - g) / (3 + 2 * (1 - g))
    c3 <- (3 - g) / (3 + 2 * (1 - g))
    EE <- 1 + A2
    VV <- 2 * (1 + B)
    EEstar <- 1 / (1 - A2)
    VVstar <- 2 * ((1 + c1 * B) / ((1 - c2 * B)^2 * (1 - c3 * B)))
    V0 <- 1 + c1 * B
    V1 <- 1 - c2 * B
    V2 <- 1 - c3 * B
    V0 <- ifelse(abs(V0) < 1e-10, 0, V0)
    rho <- (.divZero(1 - A2, V1))^2 * V0 / V2
    df2 <- 4 + 3 / (rho - 1)
    df2
}

.divZero <- function(x, y, tol = 1e-14) {
    ## ratio x/y is set to 1 if both |x| and |y| are below tol
    if (abs(x) < tol && abs(y) < tol) {
        1
    } else {
        x / y
    }
}

## FIXME: why do we go through this?
.vcov_kenward_adjusted <- function(model) {
    .vcovAdj16_internal(stats::vcov(model)$cond, .get_SigmaG(model), glmmTMB::getME(model, "X"))
}

.get_SigmaG <- function(model) {

    GGamma <- VarCorr(model)$cond
    SS <- .shgetME(model)

    ## Put covariance parameters for the random effects into a vector:
    ## TODO: It is a bit ugly to throw everything into one long vector here; a list would be more elegant
    ggamma <- list()
    for (ii in 1:(SS$n.RT)) {
        Lii <- GGamma[[ii]]
        ggamma <- c(ggamma, Lii[lower.tri(Lii, diag = TRUE)])
    }
    ggamma[[length(ggamma)+1]] <- family(model)$linkfun(.family_var_func(model)) ## Extend ggamma by the residuals variance
    n.ggamma <- length(ggamma)

    ## Find G_r:
    G <- NULL
    Zt <- Matrix::t(getME(model, "Z"))
    for (ss in 1:SS$n.RT) {
        ZZ <- .shget_Zt_group(ss, Zt, SS$Gp)
        n.lev <- SS$n.lev.by.RT2[ss] ## ; cat(sprintf("n.lev=%i\n", n.lev))
        Ig <- Matrix::sparseMatrix(1:n.lev, 1:n.lev, x = 1)
        for (rr in 1:SS$n.parm.by.RT[ss]) {
            ## This is takes care of the case where there is random regression and several matrices have to be constructed.
            ## FIXME: I am not sure this is correct if there is a random quadratic term. The '2' below looks suspicious.
            ii.jj <- .index2UpperTriEntry(rr, SS$n.comp.by.RT[ss]) ## ; cat("ii.jj:"); print(ii.jj)
            ii.jj <- unique(ii.jj)
            if (length(ii.jj) == 1) {
                EE <- Matrix::sparseMatrix(
                                  ii.jj,
                                  ii.jj,
                                  x = 1,
                                  dims = rep(SS$n.comp.by.RT[ss], 2)
                              )
            } else {
                EE <- Matrix::sparseMatrix(ii.jj, ii.jj[2:1], dims = rep(SS$n.comp.by.RT[ss], 2))
            }
            EE <- Ig %x% EE ## Kronecker product
            G <- c(G, list(t(ZZ) %*% EE %*% ZZ))
        }
    }

    ## Extend by the identity for the residual
    n.obs <- nobs(model)
    G <- c(G, list(Matrix::sparseMatrix(1:n.obs, 1:n.obs, x = 1)))

    Sigma <- ggamma[[1]] * G[[1]]
    for (ii in 2:n.ggamma) {
        Sigma <- Sigma + ggamma[[ii]] * G[[ii]]
    }

    list(Sigma = Sigma, G = G, n.ggamma = n.ggamma)
}

.index2UpperTriEntry <- function(k, N) {
    ## inverse of indexSymmat2vec
    ## result: index pair (i,j) with i>=j
    ## k: element in the vector of upper triangular elements
    ## example: N=3: k=1 -> (1,1), k=2 -> (1,2), k=3 -> (1,3), k=4 -> (2,2)
    aa <- cumsum(N:1)
    aaLow <- c(0, aa[-length(aa)])
    i <- which(aaLow < k & k <= aa)
    j <- k - N * i + N - i * (3 - i) / 2 + i
    c(i, j)
}

.vcovAdj16_internal <- function(Phi, SigmaG, X) {

    ## slow (inverse of nxn matrix)
    SigmaInv <- chol2inv(chol(Matrix::forceSymmetric(as.matrix(SigmaG$Sigma))))
    n.ggamma <- SigmaG$n.ggamma
    TT <- as.matrix(SigmaInv %*% X)
    HH <- OO <- vector("list", n.ggamma)

    for (ii in 1:n.ggamma) {
        HH[[ii]] <- as.matrix(SigmaG$G[[ii]] %*% SigmaInv)
        OO[[ii]] <- as.matrix(HH[[ii]] %*% X)
    }

    ## Finding PP, QQ
    PP <- QQ <- NULL
    for (rr in 1:n.ggamma) {
        OrTrans <- t(OO[[rr]])
        PP <- c(PP, list(Matrix::forceSymmetric(-1 * OrTrans %*% TT)))
        for (ss in rr:n.ggamma) {
            QQ <- c(QQ, list(OrTrans %*% SigmaInv %*% OO[[ss]]))
        }
    }

    PP <- as.matrix(PP)
    QQ <- as.matrix(QQ)

    Ktrace <- matrix(NA, nrow = n.ggamma, ncol = n.ggamma)
    for (rr in 1:n.ggamma) {
        HrTrans <- t(HH[[rr]])
        for (ss in rr:n.ggamma) {
            Ktrace[rr, ss] <- Ktrace[ss, rr] <- sum(HrTrans * HH[[ss]])
        }
    }

    ## Finding information matrix
    IE2 <- matrix(NA, nrow = n.ggamma, ncol = n.ggamma)

    for (ii in 1:n.ggamma) {
        Phi.P.ii <- Phi %*% PP[[ii]]
        for (jj in ii:n.ggamma) {
            www <- .indexSymmat2vec(ii, jj, n.ggamma)
            IE2[ii, jj] <- IE2[jj, ii] <- Ktrace[ii, jj] -
                2 * sum(Phi * QQ[[www]]) + sum(Phi.P.ii * (PP[[jj]] %*% Phi))
        }
    }

    eigenIE2 <- eigen(IE2, only.values = TRUE)$values
    condi <- min(abs(eigenIE2))

    WW <- if (condi > 1e-10) {
              as.matrix(Matrix::forceSymmetric(2 * solve(IE2)))
          } else {
              as.matrix(Matrix::forceSymmetric(2 * MASS::ginv(IE2)))
          }

    UU <- matrix(0, nrow = ncol(X), ncol = ncol(X))
    for (ii in 1:(n.ggamma - 1)) {
        for (jj in (ii + 1):n.ggamma) {
            www <- .indexSymmat2vec(ii, jj, n.ggamma)
            UU <- UU + WW[ii, jj] * (QQ[[www]] - PP[[ii]] %*% Phi %*% PP[[jj]])
        }
    }

    UU <- as.matrix(UU)
    UU <- UU + t(UU)
    for (ii in 1:n.ggamma) {
        www <- .indexSymmat2vec(ii, ii, n.ggamma)
        UU <- UU + WW[ii, ii] * (QQ[[www]] - PP[[ii]] %*% Phi %*% PP[[ii]])
    }

    GGAMMA <- Phi %*% UU %*% Phi
    PhiA <- Phi + 2 * GGAMMA
    attr(PhiA, "P") <- PP
    attr(PhiA, "W") <- WW
    attr(PhiA, "condi") <- condi
    PhiA
}



.indexSymmat2vec <- function(i, j, N) {
    ## S[i,j] symetric N times N matrix
    ## r the vector of upper triangular element  in row major order:
    ## r= c(S[1,1],S[1,2]...,S[1,j], S[1,N], S[2,2],...S[N,N]
    ## Result: k: index of k-th element of r
    k <- if (i <= j) {
             (i - 1) * (N - i / 2) + j
         } else {
             (j - 1) * (N - j / 2) + i
         }
}

.shgetME <- function(model) {
    ##Gp <- lme4::getME(model, "Gp")
    Gp <- unname(cumsum(c(0,sapply(model$modelInfo$reStruc$condReStruc, function(x) x$blockReps*x$blockSize))))
    n.RT <- length(Gp) - 1 ## Number of random terms (i.e. of (|)'s)
    n.lev.by.RT <- sapply(model$modelInfo$reTrms$cond$flist, nlevels)
    n.comp.by.RT <- .get.RT.dim.by.RT(model)
    n.parm.by.RT <- (n.comp.by.RT + 1) * n.comp.by.RT / 2
    n.RE.by.RT <- diff(Gp)

    n.lev.by.RT2 <- n.RE.by.RT / n.comp.by.RT ## Same as n.lev.by.RT2 ???

    list(
        Gp = Gp, ## group.index
        n.RT = n.RT, ## n.groupFac
        n.lev.by.RT = n.lev.by.RT, ## nn.groupFacLevelsNew
        n.comp.by.RT = n.comp.by.RT, ## nn.GGamma
        n.parm.by.RT = n.parm.by.RT, ## mm.GGamma
        n.RE.by.RT = n.RE.by.RT, ## ... Not returned before
        n.lev.by.RT2 = n.lev.by.RT2, ## nn.groupFacLevels
        n_rtrms = length(names(glmmTMB::ranef(model)))
    )
}

## Alternative to .get_Zt_group
.shget_Zt_group <- function(ii.group, Zt, Gp, ...) {
    zIndex.sub <- (Gp[ii.group] + 1):Gp[ii.group + 1]
    as.matrix(Zt[zIndex.sub, ])
}

.get.RT.dim.by.RT <- function(model) {
    ## output: dimension (no of columns) of covariance matrix for random term ii
    lengths(lapply(glmmTMB::ranef(model)$cond, colnames))
}

## Precompute (and cache) the pieces needed for Satterthwaite denominator-df
## calculations that depend only on the fitted model, not on the specific
## contrast(s) being tested: the (inverse) Hessian of the negative
## log-likelihood with respect to the variance/dispersion parameters
## ("kappa"), and the Jacobian of cov(beta) with respect to those same
## parameters.
##
## Both require repeated evaluation of expensive functions (`model$obj$gr()`,
## and -- via `.covbeta_kappa()` -- `TMB::sdreport()`) at perturbed parameter
## values, so:
##  (1) the result is cached on `model$obj$env`, a genuine R environment
##      (unlike the rest of `model`, which is an ordinary list and so is
##      copied rather than shared when passed around) that is shared by
##      reference across every copy of `model`; this means repeated calls to
##      `dof_satt()` on the same fitted model (e.g. from more than one
##      `summary(fit, ddf = "satterthwaite")` call, or from any future code
##      path -- such as a joint/multi-model Satterthwaite test -- that needs
##      the same per-model quantities) reuse this cache instead of redoing
##      the work; and
##  (2) `numDeriv::jacobian()` is called with `method = "simple"` (one-sided
##      differences, ~p+1 evaluations) rather than the default
##      `"Richardson"` (which redoes each finite difference at several step
##      sizes for extrapolated accuracy, at roughly 4-8x the function
##      evaluations); the resulting ddf are an approximation in any case, and
##      `method.args` tuning of the default Richardson method previously
##      found no detectable precision benefit (see the removed
##      `method.args = list(r = 6)` experiment below `.get_jac_list()`), so
##      there is little accuracy to lose by using cheaper differencing.
.satt_precompute <- function(model) {
    cache_env <- model$obj$env
    if (!is.null(cache_env$.satt_cache)) {
        return(cache_env$.satt_cache)
    }
    kappa_opt <- model$fit$par
    h_kappa <- numDeriv::jacobian(func = model$obj$gr, x = kappa_opt, method = "simple")
    ## one-sided ("simple") finite differences are not guaranteed to give a
    ## numerically symmetric matrix; eigen(symmetric=TRUE) would silently use
    ## only the lower triangle in that case, distorting the result, so
    ## symmetrize explicitly first
    h_kappa <- (h_kappa + t(h_kappa)) / 2
    eig_h_kappa <- eigen(h_kappa, symmetric = TRUE)
    ## diag(x) for a length-1 numeric x builds an x-by-x identity matrix
    ## rather than a 1x1 matrix containing x -- an easy footgun when there's
    ## only one outer/kappa parameter (e.g. a single random-intercept
    ## variance with no dispersion parameter to estimate, as for a poisson
    ## GLMM); nrow= makes this robust regardless of length(values)
    cov_varpar_kappa <- with(eig_h_kappa,
                             vectors %*% diag(1/values, nrow = length(values)) %*% t(vectors))
    jac_kappa <- .get_jac_list(.covbeta_kappa, kappa_opt, model, method = "simple")
    res <- list(cov_varpar_kappa = cov_varpar_kappa, jac_kappa = jac_kappa)
    cache_env$.satt_cache <- res
    res
}

#' @rdname dof_KR
#'
#' @export
#' @param L a contrast matrix: by default, equal to an identity matrix (i.e., ddfs are returned
#' for each fixed-effect parameter)
dof_satt <- function(model, L = diag(length(fixef(model)$cond))) {
    model_vcov <- vcov(model, full = TRUE)

    pre <- .satt_precompute(model)
    cov_varpar_kappa <- pre$cov_varpar_kappa
    jac_kappa <- pre$jac_kappa

    res <- numeric(nrow(L))
    for (i in seq_along(res)) {
        grad_kappa <- .get_gradient(jac_kappa, L[i,])
        var_Lbeta <- drop(t(L[i,]) %*% vcov(model)$cond %*% L[i,])
        v_numerator <- 2 * var_Lbeta ^ 2
        v_denominator_kappa <- sum(grad_kappa * (cov_varpar_kappa %*% grad_kappa))
        res[i] <- v_numerator/v_denominator_kappa
    }
    res
}

.covbeta_kappa <- function(kappa,md) {
  sdr <- TMB::sdreport(
    md$obj,
    par.fixed = kappa,
    getJointPrecision = TRUE
  )
  q_mat <- sdr$jointPrecision
  which_fixed <- which(rownames(q_mat) == "beta")
  q_marginal <- unname(GMRFmarginal(q_mat, which_fixed))
  solve(as.matrix(q_marginal))
}


.get_jac_list <- function(covbeta_fun, x_opt, md, ...) {
    jac_matrix <- numDeriv::jacobian(
                                func = covbeta_fun,
                                x = x_opt,
                                md=md,
                                ## does not help anything -
                                ## but seems precision is already good enough:
                                ## method.args = list(r = 6),
                                ...
                            )
    res <- list()
    for (i in seq_len(ncol(jac_matrix))) { # for each variance parameter
        jac_col <- jac_matrix[, i]
        p <- sqrt(length(jac_col))
        ## get p x p matrix
        res[[i]] <- matrix(jac_col, nrow = p, ncol = p)
    }
    res
}


.get_gradient <- function(jac, L) {
    vapply(
        jac,
        FUN = function(x) sum(L * x %*% L), # = {L' Jac L}_i
        FUN.VALUE = numeric(1L)
    )
}

## ---- joint (multi-parameter) ddf calculations ----
## used by anova.glmmTMB (ddf != "asymptotic") to compare two nested models
## via an F-ratio test rather than a likelihood ratio test.
##
## `.KR_adjust_joint` generalizes `.adjusted_ddf` (above) from a single contrast
## vector to a q-row contrast matrix L, and additionally returns the F-statistic
## and p-value for the joint test; adapted from the (unexported) `.KR_adjust`
## function in pbkrtest (which is itself model-class-agnostic, unlike the rest
## of pbkrtest's Kenward-Roger machinery)
.KR_adjust_joint <- function(adjusted_vcov, unadjusted_vcov, L, beta, betaH = 0) {
    Theta <- t(L) %*% solve(L %*% unadjusted_vcov %*% t(L), L)
    P <- attr(adjusted_vcov, "P")
    W <- attr(adjusted_vcov, "W")
    A1 <- A2 <- 0
    ThetaPhi <- Theta %*% unadjusted_vcov
    n.ggamma <- length(P)
    for (ii in 1:n.ggamma) {
        for (jj in ii:n.ggamma) {
            e <- if (ii == jj) 1 else 2
            ui <- ThetaPhi %*% P[[ii]] %*% unadjusted_vcov
            uj <- ThetaPhi %*% P[[jj]] %*% unadjusted_vcov
            A1 <- A1 + e * W[ii, jj] * (sum(diag(ui)) * sum(diag(uj)))
            A2 <- A2 + e * W[ii, jj] * sum(ui * t(uj))
        }
    }
    q <- as.numeric(Matrix::rankMatrix(L))
    B <- (A1 + 6 * A2) / (2 * q)
    g <- ((q + 1) * A1 - (q + 4) * A2) / ((q + 2) * A2)
    c1 <- g / (3 * q + 2 * (1 - g))
    c2 <- (q - g) / (3 * q + 2 * (1 - g))
    c3 <- (q + 2 - g) / (3 * q + 2 * (1 - g))
    V0 <- 1 + c1 * B
    V1 <- 1 - c2 * B
    V2 <- 1 - c3 * B
    V0 <- ifelse(abs(V0) < 1e-10, 0, V0)
    rho <- (.divZero(1 - A2 / q, V1))^2 * V0 / (q * V2)
    df2 <- 4 + (q + 2) / (q * rho - 1)
    F.scaling <- if (abs(df2 - 2) < 0.01) 1 else df2 * (1 - A2 / q) / (df2 - 2)

    betaDiff <- cbind(beta - betaH)
    Lb2 <- L %*% betaDiff
    Wald <- as.numeric(t(Lb2) %*% solve(L %*% adjusted_vcov %*% t(L), Lb2))
    Fstat <- F.scaling * (Wald / q)
    list(Fstat = Fstat, ndf = q, ddf = df2,
         p.value = stats::pf(Fstat, df1 = q, df2 = df2, lower.tail = FALSE))
}

##' Kenward-Roger F-ratio test comparing two nested \code{glmmTMB} models
##'
##' @param largeModel the model with more (conditional) fixed-effect parameters
##' @param smallModel the model with fewer fixed-effect parameters, nested in \code{largeModel}
##' @param betaH null-hypothesis value(s) for the restricted parameters (default 0)
##' @return a list with elements \code{Fstat}, \code{ndf}, \code{ddf}, \code{p.value}
##' @noRd
.joint_ddf_KR <- function(largeModel, smallModel, betaH = 0) {
    L <- as.matrix(pbkrtest::make_restriction_matrix(getME(largeModel, "X"),
                                                      getME(smallModel, "X")))
    adjusted_vcov <- .vcov_kenward_adjusted(largeModel)
    unadjusted_vcov <- stats::vcov(largeModel)$cond
    beta <- fixef(largeModel)$cond
    .KR_adjust_joint(adjusted_vcov, unadjusted_vcov, L, beta, betaH)
}

## combine per-eigenvalue Satterthwaite dfs into a single ddf for a joint
## (multi-parameter) test; adapted from the (unexported) `get_Fstat_ddf`
## function in pbkrtest
.combine_ddf_satt <- function(nu, tol = 1e-8) {
    if (length(nu) == 1) return(nu)
    if (all(abs(diff(nu)) < tol)) return(mean(nu))
    if (any(nu <= 2)) return(2)
    E <- sum(nu / (nu - 2))
    2 * E / (E - length(nu))
}

## Satterthwaite F-ratio test for an arbitrary hypothesis L %*% beta = betaH;
## split out from .joint_ddf_satt() so callers that already have a hypothesis
## matrix in hand (car::Anova()'s Type II/III per-term tests) don't need to
## construct it via a pair of nested models -- mirrors .KR_adjust_joint(),
## which already takes L directly rather than two models
##' @param model a fitted glmmTMB model
##' @param L a hypothesis matrix (ncol == number of conditional fixed-effect parameters)
##' @param betaH null-hypothesis value(s) for \code{L \%*\% beta} (default 0)
##' @param eps eigenvalue tolerance (relative to the largest eigenvalue), below which
##' a contrast direction is dropped from the test
##' @noRd
.satt_adjust_joint <- function(model, L, betaH = 0, eps = sqrt(.Machine$double.eps)) {
    beta <- fixef(model)$cond
    vcov_beta <- stats::vcov(model)$cond

    ## reuse the per-model cache from .satt_precompute() (shared with
    ## dof_satt()) instead of redoing the expensive kappa Hessian/Jacobian
    ## computation from scratch on every call -- car::Anova(..., ddf =
    ## "satterthwaite") calls this once per term, so without the cache
    ## (and its cheaper "simple"-differencing settings) this would be much
    ## slower than necessary for models with several terms
    pre <- .satt_precompute(model)
    cov_varpar_kappa <- pre$cov_varpar_kappa
    jac_kappa <- pre$jac_kappa

    vcov_Lbeta <- L %*% vcov_beta %*% t(L)
    eig <- eigen(vcov_Lbeta)
    d <- eig$values
    tol <- max(eps * d[1], 0)
    qq <- sum(d > tol)
    PtL <- crossprod(eig$vectors, L)[seq_len(qq), , drop = FALSE]

    betaDiff <- beta - betaH
    t2 <- drop(PtL %*% betaDiff)^2 / d[seq_len(qq)]
    Fstat <- sum(t2) / qq

    nu_m <- vapply(seq_len(qq), function(m) {
        grad_kappa <- .get_gradient(jac_kappa, PtL[m, ])
        2 * d[m]^2 / sum(grad_kappa * (cov_varpar_kappa %*% grad_kappa))
    }, numeric(1))
    ddf <- .combine_ddf_satt(nu_m)

    list(Fstat = Fstat, ndf = qq, ddf = ddf,
         p.value = stats::pf(Fstat, df1 = qq, df2 = ddf, lower.tail = FALSE))
}

##' Satterthwaite F-ratio test comparing two nested \code{glmmTMB} models
##' @inheritParams .joint_ddf_KR
##' @param eps eigenvalue tolerance (relative to the largest eigenvalue), below which
##' a contrast direction is dropped from the test
##' @noRd
.joint_ddf_satt <- function(largeModel, smallModel, betaH = 0, eps = sqrt(.Machine$double.eps)) {
    L <- as.matrix(pbkrtest::make_restriction_matrix(getME(largeModel, "X"),
                                                      getME(smallModel, "X")))
    .satt_adjust_joint(largeModel, L, betaH = betaH, eps = eps)
}

## classical (exact, for Gaussian fixed-effect-only fits) multi-parameter Wald
## F-test, used in place of `.joint_ddf_KR`/`.joint_ddf_satt` when neither model
## has random effects: there is then no variance-component uncertainty for
## Kenward-Roger/Satterthwaite to correct for, and the denominator df is simply
## the residual df (nobs - npar) of the fuller model
##' @inheritParams .joint_ddf_KR
##' @noRd
## classical Wald F-test for L %*% beta = betaH at a caller-supplied ddf
## (no Kenward-Roger/Satterthwaite correction); factored out so both
## .joint_ddf_none() (model-comparison form) and car::Anova()'s per-term
## no-random-effects fallback can share it
##' @noRd
.wald_joint_test <- function(unadjusted_vcov, L, beta, ddf, betaH = 0) {
    q <- as.numeric(Matrix::rankMatrix(L))
    Lb2 <- L %*% cbind(beta - betaH)
    Wald <- as.numeric(t(Lb2) %*% solve(L %*% unadjusted_vcov %*% t(L), Lb2))
    Fstat <- Wald / q
    list(Fstat = Fstat, ndf = q, ddf = ddf,
         p.value = stats::pf(Fstat, df1 = q, df2 = ddf, lower.tail = FALSE))
}

.joint_ddf_none <- function(largeModel, smallModel, betaH = 0) {
    L <- as.matrix(pbkrtest::make_restriction_matrix(getME(largeModel, "X"),
                                                      getME(smallModel, "X")))
    .wald_joint_test(stats::vcov(largeModel)$cond, L, fixef(largeModel)$cond,
                      ddf = stats::df.residual(largeModel), betaH = betaH)
}

## the Kenward-Roger correction is derived from REML variance-component
## estimates; there is no valid correction to compute for an ML fit, so
## (unlike the softer warnings elsewhere) this is a hard error rather than
## a warning that leaves downstream code to silently proceed anyway.
## Kept in one place (rather than inlined separately in check_ddf() and
## emm_basis.glmmTMB()) so the error text can't drift out of sync between
## summary()/anova() and emmeans().
##' @noRd
.check_KR_reml <- function(object) {
    if (!isREML(object)) {
        stop("ddf='kenward-roger' requires a REML fit (fit with REML=TRUE)", call. = FALSE)
    }
    invisible(NULL)
}

## shared wording (used by check_ddf() and emm_basis.glmmTMB()) for the
## soft warning issued when Kenward-Roger/Satterthwaite are used on a
## non-Gaussian family, where their behavior is not well studied
##' @noRd
.warn_ddf_glmm <- function(ddf) {
    warning(sprintf(
        "performance (and theoretical justification) of ddf='%s' for GLMMs is poorly understood",
        ddf))
}

##' check whether a requested ddf choice is valid/sensible for a given model,
##' issuing warnings/messages (or an error, if the choice is unsupported for
##' this model's family) as needed; shared by \code{summary.glmmTMB} and
##' \code{anova.glmmTMB}
##' @param object a fitted glmmTMB model
##' @param ddf ddf choice, as in \code{summary.glmmTMB}
##' @noRd
check_ddf <- function(object, ddf) {
    if (ddf == "asymptotic") return(invisible(NULL))
    if (!hasRandom(object)) {
        ## the residual-df Wald F-test fallback used here (see
        ## .joint_ddf_none()/.wald_joint_test(), and the "!hasRandom"
        ## branches of summary.glmmTMB()/Anova.glmmTMB()) doesn't touch any
        ## family-specific machinery, so it works the same regardless of
        ## family; no need for the checks below
        message(sprintf(
            "no random effects in model: kenward-roger/satterthwaite corrections are not meaningful; using residual df (nobs - npar) for ddf='%s' instead",
            ddf))
        return(invisible(NULL))
    }
    if (ddf == "kenward-roger") {
        .check_KR_reml(object)
        ## the Kenward-Roger variance-component machinery (.get_SigmaG(),
        ## via .family_var_func()) only has cases for a handful of families
        ## with an estimated dispersion parameter (gaussian, nbinom1/2/12,
        ## Gamma); families with no dispersion parameter at all (binomial,
        ## poisson, ...) are never supported and fail with an opaque error
        ## deep inside that machinery, so reject them here with a clear
        ## message instead. Satterthwaite has no such restriction: it works
        ## directly from the TMB joint precision matrix and is family-agnostic
        ## (confirmed to work for e.g. a poisson GLMM)
        if (!usesDispersion(family(object)$family)) {
            stop(sprintf(
                "ddf='kenward-roger' is not supported for family '%s' (no estimated dispersion parameter); use ddf='satterthwaite' or ddf='asymptotic' instead",
                family(object)$family), call. = FALSE)
        }
    }
    ## applies to both kenward-roger and satterthwaite: fires after (not
    ## instead of) the checks above, so an ML fit or an unsupported family
    ## gets only that more specific error, not a spurious extra warning for
    ## a test that's about to fail anyway
    if (family(object)$family != "gaussian") {
        .warn_ddf_glmm(ddf)
    }
    if (ddf == "kenward-roger" && (!trivialDisp(object) || !noZI(object))) {
        message("ddf='kenward-roger' ignored except for conditional-distribution parameters")
    }
    invisible(NULL)
}
