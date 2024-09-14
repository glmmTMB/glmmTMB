## https://github.com/glmmTMB/glmmTMB/issues/1059#issuecomment-2187660510
library(nloptr)
library(glmmTMB)

nl_algs <- (nloptr.get.default.options()[1,"possible_values"]
    |> strsplit(",")
    |> unlist()
    |> trimws()
    |> grep(pattern="NLOPT_LD", value = TRUE)
    |> grep(pattern="AUGLAG", value = TRUE, invert = TRUE)
)

nloptwrap1 <- function(start, objective, gradient, control,
                         algorithm = "NLOPT_LD_LBFGS") {
    ## what are nlminb default tolerances?
    my_obj <- function(x) objective(x) ##
    my_grad <- function(x) gradient(x) ##
    fit <- nloptr(x0 = start, eval_f = my_obj, eval_grad_f = my_grad,
                  opts = list(algorithm = algorithm, xtol_rel = 1e-8))
    return(fit)
}

optim_algs <- c("BFGS", "CG", "L-BFGS-B")
allFit.methods <- data.frame(
    optimizer = rep(c("nlminb", "optim", "nloptwrap1"),
                    times = c(1, length(optim_algs), length(nl_algs))),
    alg = c("default", optim_algs, nl_algs))

## FIXME: don't write over existing control ('smart_refit' branch has better
##        storage of control info)
## FIXME: don't refit with original alg?
## FIXME: allow for setting tolerances
## FIXME: factory etc. to store warnings
## FIXME: summary/tidy methods?
## FIXME: exclude terrible/sensitive optimizers
## FIXME: expand to use optimx methods?
## FIXME: we should only re-run the optimization, not update()! -> finalizeTMB()
##        if necessary do the modular steps again, ideally re-use info from fitted model
##        optional whether we should start from best fit or from default/specified starting values
##        optional whether to include original alg
allFit <- function(fit, methods = allFit.methods) {
    res <- vector("list", length = nrow(methods))
    names(res) <- paste(methods$optimizer,
                        methods$alg,
                        sep = ".")
    for (i in seq(nrow(methods))) {
        optimizer <- methods$optimizer[i]
        alg <- methods$alg[i]
        cat("** ", i, optimizer, alg, "\n")
        ctrl <- switch(optimizer,
            nlminb = glmmTMBControl(optimizer=nlminb),
            optim = glmmTMBControl(optimizer=optim,
                                   optArgs = list(method = alg)),
            nloptwrap1 = glmmTMBControl(optimizer=nloptwrap1,
                                        optArgs = list(algorithm = alg)))
        time <- system.time(curfit <- try(update(fit, control = ctrl)))
        attr(curfit, "system.time") <- time
        res[[i]] <- curfit
    }
    res
}

if (FALSE) {
    ## minimal example for testing allFit
    x0 <- glmmTMB(mpg ~ cyl, data = mtcars)
    x1 <- update(x0, control = glmmTMBControl(optimizer = nloptwrap1))
    x2 <- update(x0, control = glmmTMBControl(optimizer = nloptwrap2),
                 optArgs = list(algorithm = nl_algs[1]))

    ## not identical (don't want them to be)
    fixef(x1)$cond-fixef(x2)$cond
}
