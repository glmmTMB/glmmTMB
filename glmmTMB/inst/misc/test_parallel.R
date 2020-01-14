## remotes::install_local(".")

parallel_test <- function(nt=NA, N=1e5, groups=200, seed=1) {
    require(glmmTMB)
    nt_orig <- nt
    if (is.na(nt)) {
        nt <- parallel::detectCores()
    } else {
        nt <- min(parallel::detectCores(),nt)
    }
    if (!is.null(seed)) set.seed(seed)
    xdata <- rnorm(N, 1, 2)
    data_use <- data.frame(obs = 1:N)
    data_use <- within(data_use,
    {
        
        group_var <- rep(seq(groups), times = nrow(data_use) / groups)
        group_intercept <- rnorm(groups, 0, 0.1)[group_var]
        xdata <- xdata
        ydata <- 0.3 + group_intercept + 0.5*xdata + rnorm(N, 0, 0.25)
    })
    (t_serial <- system.time(
         model3 <- glmmTMB(formula = ydata ~ 1 + xdata + (1 | group_var),
                           data = data_use,
                           control = glmmTMBControl(parallel = 1))
     )
    )
    (t_parallel <- system.time(
         update(model3,  control = glmmTMBControl(parallel = nt))
     )
    )
    ret <- list(t_serial=t_serial,
                t_parallel=t_parallel,
                pars=list(nt=nt, N=N, groups=groups, seed=seed),
                s_info=sessionInfo(),
                p_info=help(package="glmmTMB"))
    class(ret) <- "partest"
    return(ret)
}

print.partest <- function(x, ...) {
    s <- x$t_serial[["elapsed"]]
    p <- x$t_parallel[["elapsed"]]
    cat(sprintf(
       "elapsed time for N=%1.1g: serial=%1.1f, parallel=%1.1f (%d threads)\n",
        x$pars$N,s,p,x$pars$nt))
    cat(sprintf("ratio=%1.1f\n",s/p))
    cat("platform: ",x$s_info$platform,"\n")
}

if (FALSE) {
    ## test 
    (p <- parallel_test(N=1e4))
}
