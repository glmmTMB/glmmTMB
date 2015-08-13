## Copy-paste from cpp-source:
.valid_family <- c(poisson_family           = 0,
                   binomial_family          = 1,
                   negative_binomial_family = 2,
                   Gamma_family             = 3,
                   beta_family              = 4,
                   gaussian_family          = 5,
                   truncated_poisson_family = 6,
                   trunc_NB_family          = 7,
                   logistic_family          = 8,
                   betabinomial_family      = 9 )
.valid_link <- c(log_link                 = 0,
                 logit_link               = 1,
                 probit_link              = 2,
                 inverse_link             = 3,
                 cloglog_link             = 4,
                 identity_link            = 5 )
names(.valid_family) <- sub("_family","",names(.valid_family))
names(.valid_link) <- sub("_link","",names(.valid_link))

## TODO: Steamline construction of design matrices
## - Can we use lme4 utilities ?
glmmTMB <- function (
    obsname,   ## Temporary solution !
    formulaX,  ## Temporary solution !
    formulaZ,  ## Temporary solution !
    formulaZI, ## Temporary solution !
    data = NULL,
    family = "gaussian",
    link = "identity", ...)
{
    family <- match.arg(family, names(.valid_family))
    link <- match.arg(link, names(.valid_link))
    ## Get X, Z and ZI (TODO: Sparse ? )
    X <- model.matrix(formulaX, data)
    Z <- model.matrix(formulaZ, data)
    ZI <- model.matrix(formulaZI, data)
    ## Get observation vector
    yobs <- get(obsname, data)
    data.tmb <- list(
        X = X, Z = Z, ZI = ZI,
        yobs = yobs,
        family = .valid_family[family],
        link = .valid_link[link]
        )
    parameters <- list(
        beta     = rep(0, ncol(X)) ,
        b        = rep(0, ncol(Z)) ,
        betazi   = rep(0, ncol(ZI)),
        logsigma = rep(0, ncol(X)) ,
        logalpha = 0
        )
    obj <- MakeADFun(data.tmb,
                     parameters,
                     random = "b",
                     profile = NULL, ## TODO: Optionally "beta"
                     silent = FALSE, ## TODO: set to TRUE
                     DLL="glmmTMB")
    obj ## For now give the object without optimizing
}
