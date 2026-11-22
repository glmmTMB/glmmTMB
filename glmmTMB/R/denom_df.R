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
#'
#' Both approximations are computed for the parameters that are actually
#' estimated, and reported for every coefficient in \code{fixef(model)$cond}:
#' coefficients tied to each other via \code{map} share one value, while
#' coefficients fixed via \code{map} (known constants with zero variance) and
#' coefficients dropped from a rank-deficient model (not estimable, with
#' \code{NA} estimates) get \code{NA}. Likewise, a contrast (row of \code{L})
#' gets \code{NA} if it puts weight on a dropped coefficient or if it involves
#' only fixed coefficients.
#' @param model a fitted \code{glmmTMB} object
#' @export
## avoid conflict with insight::dof_kenward ...
dof_KR <- function(model) {
    sp <- .beta_spaces(model, "cond")
    Phi <- .Phi_est(model)
    dof <- rep(NA_real_, sp$p_nom)
    names(dof) <- sp$names_nom
    krvcov <- Phi
    if (sp$p_est > 0 && .Phi_ok(Phi, "Kenward-Roger")) {
        krvcov <- .vcov_kenward_adjusted(model, sp, Phi)
        ## one unit contrast per coefficient, carried over to the estimated
        ## space (dropped and map-fixed coefficients keep NA)
        Lest <- .lift_contrasts_est(diag(nrow = sp$p_nom), sp)
        for (i in which(Lest$ok)) dof[i] <- .adjusted_ddf(krvcov, Lest$L[i, ], Phi)
    }
    ## reported in the nominal space like vcov(include_nonest = TRUE): NA
    ## rows/columns for dropped, zero rows/columns for map-fixed coefficients
    krvcov <- .lift_vcov_nominal(krvcov, sp)
    attr(dof, "vcov") <- krvcov
    attr(dof, "se") <- abs(sqrt(diag(krvcov)))
    dof
}

## Covariance matrix of the *estimated* conditional fixed effects (see
## .beta_spaces()): what the Kenward-Roger and Satterthwaite machinery
## works with. A 0 x 0 matrix if no coefficient is estimated at all (every
## coefficient fixed via 'map'), in which case vcov() returns NULL
.Phi_est <- function(model, component = "cond") {
    Phi <- stats::vcov(model, include_nonest = FALSE)[[component]]
    if (is.null(Phi)) return(matrix(numeric(0), 0, 0))
    as.matrix(Phi)
}

## sdreport() can return an NA-filled covariance matrix for a boundary fit
## (random-effects variance -> 0) even with pdHess = TRUE; warn once and let
## the callers return NA rather than fail inside eigen()
.Phi_ok <- function(Phi, what) {
    if (!anyNA(Phi)) return(TRUE)
    warning(sprintf(
        "covariance matrix of the fixed effects contains NA values (see diagnose()); %s degrees of freedom set to NA",
        what), call. = FALSE)
    FALSE
}

## Carry contrasts (one per row of L) over to the estimated space. L may be
## given on the nominal coefficients (as summary()'s identity default), on
## the X-space (non-NA) coefficients (as emmeans supplies them) or already
## on the estimated parameters; the three can only coincide in width when
## they coincide in meaning. Returns the estimated-space matrix and a
## per-row flag that is FALSE for contrasts without a df: those with weight
## on a rank-dropped (not estimable) coefficient, and those that vanish in
## the estimated space (involving only map-fixed coefficients, so that the
## contrast is a known constant with zero variance)
.lift_contrasts_est <- function(L, sp) {
    if (!is.matrix(L)) L <- matrix(L, nrow = 1)
    ok <- rep(TRUE, nrow(L))
    if (ncol(L) == sp$p_nom && sp$p_X < sp$p_nom) {
        ok <- rowSums(L[, !sp$keep, drop = FALSE] != 0) == 0
        L <- L[, sp$keep, drop = FALSE]
    }
    if (ncol(L) == sp$p_X && sp$mapped) {
        L <- L %*% sp$A
    } else if (ncol(L) != sp$p_est) {
        stop(sprintf("contrast matrix has %d columns; expected %s (one per fixed-effect coefficient)",
                     ncol(L), paste(unique(c(sp$p_nom, sp$p_X, sp$p_est)), collapse = " or ")))
    }
    ok <- ok & rowSums(L != 0) > 0
    list(L = L, ok = ok)
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

## Kenward-Roger adjusted covariance matrix, in the estimated space (see
## .beta_spaces()): X is restricted to the estimated parameters, X %*% A,
## so that it conforms with the covariance matrix of those parameters.
## Carries the "P" and "W" attributes needed by .adjusted_ddf()
.vcov_kenward_adjusted <- function(model, sp = .beta_spaces(model, "cond"), Phi = .Phi_est(model)) {
    X <- glmmTMB::getME(model, "X")
    if (sp$mapped) X <- X %*% sp$A
    .vcovAdj16_internal(Phi, .get_SigmaG(model), X)
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
#' for each fixed-effect parameter). Columns may correspond to all coefficients in
#' \code{fixef(model)$cond}, to the non-\code{NA} (estimable) ones, or to the parameters
#' actually estimated (after \code{map})
dof_satt <- function(model, L = diag(length(fixef(model)$cond))) {
    if (!is.matrix(L)) L <- matrix(L, nrow = 1)
    sp <- .beta_spaces(model, "cond")
    Lest <- .lift_contrasts_est(L, sp)
    res <- rep(NA_real_, nrow(L))
    Phi <- .Phi_est(model)
    if (sp$p_est == 0 || !any(Lest$ok) || !.Phi_ok(Phi, "Satterthwaite")) return(res)

    pre <- .satt_precompute(model)
    cov_varpar_kappa <- pre$cov_varpar_kappa
    jac_kappa <- pre$jac_kappa

    for (i in which(Lest$ok)) {
        l <- Lest$L[i, ]
        grad_kappa <- .get_gradient(jac_kappa, l)
        var_Lbeta <- drop(t(l) %*% Phi %*% l)
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
## via an F-ratio test rather than a likelihood ratio test, and by
## car::Anova() for its per-term Type II/III F tests. The hypothesis
## L %*% beta = betaH comes with L on the X-space coefficients (what
## car::Anova() and pbkrtest::make_restriction_matrix() supply) and is
## carried over to the estimated space by .hypothesis_est() first.

## L_est = L %*% A (see .beta_spaces()). Rows that vanish there involve only
## map-fixed coefficients (known constants): untestable, dropped, as
## car::Anova() does for its Wald tests. Rows that become linearly dependent
## (e.g. two tied coefficients tested separately) are reduced to a row
## basis, which the F-test algebra requires; this is only legitimate if the
## dropped rows say the same thing as the kept ones once the fixed values
## are substituted, otherwise the hypothesis contradicts itself (e.g. a
## "nested" model that fixes a coefficient to a nonzero value) and gets NA.
## The estimate Lb = L %*% beta_X - betaH is formed on the X-space
## coefficients so that the fixed values enter it.
.hypothesis_est <- function(L, beta_X, betaH = 0, sp) {
    L <- as.matrix(L)
    if (ncol(L) != sp$p_X) {
        stop(sprintf("hypothesis matrix has %d columns; expected %d (one per estimable coefficient)",
                     ncol(L), sp$p_X))
    }
    Lb <- L %*% cbind(beta_X - betaH)
    L_est <- if (sp$mapped) L %*% sp$A else L
    ok <- rowSums(L_est != 0) > 0
    L_est <- L_est[ok, , drop = FALSE]
    Lb <- Lb[ok, , drop = FALSE]
    if (nrow(L_est) > 1 && (qrL <- qr(t(L_est)))$rank < nrow(L_est)) {
        if (max(abs(qr.resid(qr(L_est), Lb))) > 1e-8 * max(1, abs(Lb))) {
            warning("hypothesis rows contradict each other once the coefficients fixed via 'map' ",
                    "are substituted; test set to NA", call. = FALSE)
            return(list(L = L_est[0, , drop = FALSE], Lb = Lb[0, , drop = FALSE], q = 0L))
        }
        rows <- sort(qrL$pivot[seq_len(qrL$rank)])
        L_est <- L_est[rows, , drop = FALSE]
        Lb <- Lb[rows, , drop = FALSE]
    }
    list(L = L_est, Lb = Lb, q = nrow(L_est))
}

## `.KR_adjust_joint` generalizes `.adjusted_ddf` (above) from a single contrast
## vector to a q-row contrast matrix L, and additionally returns the F-statistic
## and p-value for the joint test; adapted from the (unexported) `.KR_adjust`
## function in pbkrtest (which is itself model-class-agnostic, unlike the rest
## of pbkrtest's Kenward-Roger machinery). L and both covariance matrices
## live in the estimated space (see .hypothesis_est()); Lb is the estimate
## of the hypothesis, L %*% beta - betaH
.KR_adjust_joint <- function(adjusted_vcov, unadjusted_vcov, L, Lb) {
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
    q <- nrow(L)
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

    Wald <- as.numeric(t(Lb) %*% solve(L %*% adjusted_vcov %*% t(L), Lb))
    Fstat <- F.scaling * (Wald / q)
    list(Fstat = Fstat, ndf = q, ddf = df2,
         p.value = stats::pf(Fstat, df1 = q, df2 = df2, lower.tail = FALSE))
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

## Satterthwaite F-ratio test for a hypothesis given in the estimated space
## (L, and its estimate Lb = L %*% beta - betaH, from .hypothesis_est());
## mirrors .KR_adjust_joint(), which likewise takes L directly rather than
## two models
##' @param model a fitted glmmTMB model
##' @param L hypothesis matrix in the estimated space
##' @param Lb estimate of the hypothesis, \code{L \%*\% beta - betaH}
##' @param vcov_beta covariance matrix of the estimated coefficients
##' @param eps eigenvalue tolerance (relative to the largest eigenvalue), below which
##' a contrast direction is dropped from the test
##' @noRd
.satt_adjust_joint <- function(model, L, Lb, vcov_beta, eps = sqrt(.Machine$double.eps)) {
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

    t2 <- drop(crossprod(eig$vectors, Lb))[seq_len(qq)]^2 / d[seq_len(qq)]
    Fstat <- sum(t2) / qq

    nu_m <- vapply(seq_len(qq), function(m) {
        grad_kappa <- .get_gradient(jac_kappa, PtL[m, ])
        2 * d[m]^2 / sum(grad_kappa * (cov_varpar_kappa %*% grad_kappa))
    }, numeric(1))
    ddf <- .combine_ddf_satt(nu_m)

    list(Fstat = Fstat, ndf = qq, ddf = ddf,
         p.value = stats::pf(Fstat, df1 = qq, df2 = ddf, lower.tail = FALSE))
}

## classical Wald F-test for L %*% beta = betaH at a caller-supplied ddf (no
## Kenward-Roger/Satterthwaite correction), used when the model has no random
## effects: there is then no variance-component uncertainty to correct for,
## and the denominator df is simply the residual df (nobs - npar). L, Lb and
## the covariance matrix live in the estimated space, as above
##' @noRd
.wald_joint_test <- function(unadjusted_vcov, L, Lb, ddf) {
    q <- nrow(L)
    Wald <- as.numeric(t(Lb) %*% solve(L %*% unadjusted_vcov %*% t(L), Lb))
    Fstat <- Wald / q
    list(Fstat = Fstat, ndf = q, ddf = ddf,
         p.value = stats::pf(Fstat, df1 = q, df2 = ddf, lower.tail = FALSE))
}

## Common driver for the joint F tests (ddf = "kenward-roger",
## "satterthwaite" or "none" = Wald F with residual df, which is also what
## a model without random effects gets): L on the X-space coefficients. A
## hypothesis with nothing testable left, or a model with no estimated
## coefficient or an NA-filled covariance matrix, gives an NA row rather
## than an error. 'info' (from .joint_test_setup()) lets callers reuse the
## model-level pieces, notably the Kenward-Roger adjusted vcov, across
## several hypotheses
.joint_test <- function(model, L, ddf, betaH = 0, info = .joint_test_setup(model, ddf)) {
    empty <- list(Fstat = NA_real_, ndf = 0, ddf = NA_real_, p.value = NA_real_)
    sp <- info$sp
    if (sp$p_est == 0) return(empty)
    hyp <- .hypothesis_est(L, fixef(model)$cond[sp$keep], betaH, sp)
    if (hyp$q == 0) return(empty)
    if (!info$Phi_ok) {
        empty$ndf <- hyp$q
        return(empty)
    }
    switch(info$ddf,
           "kenward-roger" = .KR_adjust_joint(info$adjusted_vcov, info$Phi, hyp$L, hyp$Lb),
           "satterthwaite" = .satt_adjust_joint(model, hyp$L, hyp$Lb, info$Phi),
           "none" = .wald_joint_test(info$Phi, hyp$L, hyp$Lb, ddf = stats::df.residual(model)),
           stop(sprintf("unknown ddf specification '%s'", info$ddf)))
}

.joint_test_setup <- function(model, ddf) {
    if (!hasRandom(model)) ddf <- "none"
    sp <- .beta_spaces(model, "cond")
    Phi <- .Phi_est(model)
    Phi_ok <- sp$p_est > 0 && .Phi_ok(Phi, "F-test denominator")
    list(ddf = ddf, sp = sp, Phi = Phi, Phi_ok = Phi_ok,
         adjusted_vcov = if (ddf == "kenward-roger" && Phi_ok) .vcov_kenward_adjusted(model, sp, Phi) else NULL)
}

##' F-ratio tests comparing two nested \code{glmmTMB} models
##'
##' @param largeModel the model with more (conditional) fixed-effect parameters
##' @param smallModel the model with fewer fixed-effect parameters, nested in \code{largeModel}
##' @param betaH null-hypothesis value(s) for the restricted parameters (default 0)
##' @return a list with elements \code{Fstat}, \code{ndf}, \code{ddf}, \code{p.value}
##' @noRd
.joint_ddf_models <- function(largeModel, smallModel, ddf, betaH = 0) {
    L <- as.matrix(pbkrtest::make_restriction_matrix(getME(largeModel, "X"),
                                                      getME(smallModel, "X")))
    .joint_test(largeModel, L, ddf, betaH = betaH)
}
.joint_ddf_KR <- function(largeModel, smallModel, betaH = 0) {
    .joint_ddf_models(largeModel, smallModel, "kenward-roger", betaH)
}
.joint_ddf_satt <- function(largeModel, smallModel, betaH = 0) {
    .joint_ddf_models(largeModel, smallModel, "satterthwaite", betaH)
}
## used in place of `.joint_ddf_KR`/`.joint_ddf_satt` when neither model
## has random effects
.joint_ddf_none <- function(largeModel, smallModel, betaH = 0) {
    .joint_ddf_models(largeModel, smallModel, "none", betaH)
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
