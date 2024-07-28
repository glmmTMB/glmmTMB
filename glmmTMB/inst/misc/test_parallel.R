## remotes::install_local(".")
cc <- commandArgs(trailingOnly = TRUE)
do_test <- FALSE
if (length(cc) > 0) do_test <- as.logical(cc[1])

##' @param nt number of threads (NA -> use detectCores() to use all cores)
##' @param N number of obs in sim data
##' @param group number of groups in sim data (groups/N must be integer)
##' @param seed random-number seed
parallel_test <- function(nt=NA, N=1e5, groups=200, seed=1) {
    require(glmmTMB)
    nt_orig <- nt
    if (is.na(nt)) {
        nt <- parallel::detectCores()
    } else {
        nt <- min(parallel::detectCores(), nt)
    }
    if (!is.null(seed)) set.seed(seed)
    xdata <- rnorm(N, 1, 2)
    data_use <- data.frame(obs = 1:N)
    data_use <- within(data_use, {
        group_var <- rep(seq(groups), times = nrow(data_use) / groups)
        group_intercept <- rnorm(groups, 0, 0.1)[group_var]
        xdata <- xdata
        ydata <- 0.3 + group_intercept + 0.5*xdata + rnorm(N, 0, 0.25)
    })
    t_serial <- system.time(
        model3 <- glmmTMB(formula = ydata ~ 1 + xdata + (1 | group_var),
                          data = data_use,
                          control = glmmTMBControl(parallel = 1))
    )
    t_parallel_noauto <- system.time(
        update(model3,  control = glmmTMBControl(parallel = list(n = nt, autopar = FALSE)))
    )
    t_parallel_auto <- system.time(
        update(model3,  control = glmmTMBControl(parallel = list(n = nt, autopar = TRUE)))
    )

    ret <- list(t_serial=t_serial,
                t_parallel_noauto=t_parallel_noauto,
                t_parallel_auto  =t_parallel_auto,
                pars=list(nt=nt, N=N, groups=groups, seed=seed),
                s_info=sessionInfo(),
                p_info=help(package="glmmTMB"))
    class(ret) <- "partest"
    return(ret)
}

print.partest <- function(x, ...) {
    s <- x$t_serial[["elapsed"]]
    pn <- x$t_parallel_noauto[["elapsed"]]
    pa <- x$t_parallel_auto[["elapsed"]]
    cat(sprintf(
        "elapsed time for N=%1.1g: serial=%1.1f, parallel_noauto=%1.1f, parallel_auto=%1.1f (%d threads)\n",
        x$pars$N,s,pn,pa,x$pars$nt))
    cat(sprintf("ratios=%1.1f (noauto), %1.1f (auto), %1.1f (noauto vs auto)\n",s/pn, s/pa, pn/pa))
    cat("platform:", x$s_info$platform, "\n")
}

## test
options(glmmTMB_openmp_debug = TRUE)
print(p <- parallel_test(nt = 4, N=1e4, groups = 250))

