#' compute Kenward-Roger degrees of freedom
#'
#' This function uses an adaptation of the machinery from the \code{pbkrtest} package
#' to compute the Kenward-Roger approximation of the 'denominator degrees of freedom' for
#' each fixed-effect coefficient in the conditional model.
#' @return a named vector of ddf for each conditional fixed-effect parameter, with attributes 'vcov'
#' (Kenward-Roger adjusted covariance matrix) and 'se' (the corresponding standard errors)
#' @details this function \emph{should not be used} for models fitted with ML rather than REML;
#' the theory is only well understood, and the model is only tested, for LMMs (\code{family = "gaussian"}).
#' Use at your own risk for GLMMs!
#' @param model a fitted \code{glmmTMB} object
#' @export
## FIXME: possible confusion/conflict with insight::dof_kenward ...
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
