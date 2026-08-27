remotes::install_github("mrc-ide/memprof")
library(memprof)
## could also try peakRAM, but doesn't seem to be reliable at
## capturing TMB memory use
library(reformulas)
library(glmmTMB)
library(dplyr)

simfun <- function(nobs=1e4, ngrp=1e2, ntimes=1e2,
                    reterm = "ou(0+ftime|g)", seed = NULL,
                    pars = list(theta=c(1,1))) {
    if (!is.null(seed)) set.seed(seed)
    dd <- data.frame(g = sample(ngrp, size = nobs, replace = TRUE),
                     time = sample(ntimes, size = nobs, replace = TRUE))
    dd$ftime <- numFactor(dd$time)
    dd$time <- factor(dd$time)
    form <- reformulate(c("1",reterm))
    dd$y <- simulate_new(form, newdata = dd, newparams = pars,
                         family = gaussian)[[1]]
    dd
}

fitfun <- function(data, reterm = "ou(0+ftime|g)") {
    form <- reformulate(c("1",reterm), response = "y")
    glmmTMB(form, data = data, family = gaussian)
}

set.seed(101)
dd <- simfun()
res <- with_monitor(fitfun(dd))
plot(rss~time, data = res$memory_use)
## max mem use, in GB
max(res$memory_use$rss/1e9) ## 0.818

dd <- simfun(nobs = 1e5, ntimes = 1e3, seed = 101)
res <- with_monitor(fitfun(dd))
max(res$memory_use$rss/1e9) ## 2.84

res_ar1 <- with_monitor(fitfun(dd, "ar1(0+time|g)"))
max(res_ar1$memory_use$rss/1e9)  ## 3.67


###
nt <- 2e4
nn <- 3e5

set.seed(101)
## counts with at least 2 times per individual
nt_id <- rpois(nt, lambda = mean(nn/nt)-2) + 2
nn2 <-sum(nt_id)  ## close to nn

tvec <- lapply(nt_id, seq.int) |> unlist() |> factor()
idvec <- rep(seq_along(nt_id), times = nt_id)
dd <- (replicate(6, rnorm(nn2), simplify = FALSE)
    |> as.data.frame()
    |> setNames(c("y", paste0("x", 1:5)))
)
dd <- data.frame(dd, tvec, id = idvec)

res <- with_monitor(
    glmmTMB(y ~ x1 + x2 + x3 + x4 + x5 + ar1(0 + tvec|id) + (1|id),
            dispformula = ~ 0,
            data = dd, verbose = TRUE)
)
