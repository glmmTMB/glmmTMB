source("misc/leverage_funs.R")
glmmTMB_lev_dir <- "glmmTMB_lev"

## reinstalling ...
if (FALSE) {

    remotes::install_github("mrc-ide/memprof")
    
    pkgdir <- "~/R/pkgs/glmmTMB/glmmTMB"
    git_status <- system(sprintf("git -C %s status", pkgdir), intern = TRUE)
    if (git_status[1] != "On branch formula_env") stop("need formula_env branch")
    install_glmmTMB(pkgdir, glmmTMB_lev_dir)
}

library(memprof)
library(glmmTMB, lib = glmmTMB_lev_dir)
library(peakRAM)
library(RTMB)
library(Matrix)
library(tidyverse)
do_simple_tests <- FALSE
do_slow_tests <- FALSE

if (do_simple_tests) {
    ## basic tests
    peakRAM_testfun(5, 5, 2)  ## 7 Mb
    peakRAM_testfun(5, 5, 3)  ## 9 Mb
    peakRAM_testfun(5, 5, 3, include_ttt = TRUE)  ## 20 Mb
}

testframe <- expand.grid(nsubj = seq(5, 40, by = 5),
                         ntax = seq(25, 200, by = 25))
write.csv(testframe, "misc/mem_testvals.csv")

testframe2 <- expand.grid(nsubj = seq(5, 40, by = 5),
                         ntax = seq(25, 200, by = 25),
                         d = 1:3,
                         include_ttt = c(TRUE, FALSE)) |>
    mutate(ntheta = ntheta(nsubj, ntax, d, include_ttt),
           nobs = ntax*nsubj) |>
    ## order by expected mem usg, to try to get in as many as possible
    ##  before out-of-memory crash
    arrange(ntheta*nobs)
write.csv(testframe2, "mem_testvals2.csv")

## first pass
if (FALSE) {
    res <- list()
    for (i in seq(nrow(testframe))) {
        nsubj <- testframe$nsubj[i]
        ntax <- testframe$ntax[i]
        cat(i, nsubj, ntax, "\n")
        res[[i]] <- with_monitor(peakRAM_testfun(nsubj, ntax))
        saveRDS(res, "memusg.rds")
    }
}

do_test <- function(i) {
    nsubj <- testframe2$nsubj[i]
    ntax <-  testframe2$ntax[i]
    d <- testframe2$d[i]
    include_ttt <- testframe2$include_ttt[i]
    ## with(testframe2[i,], {
    cat(i, nsubj, ntax, d, include_ttt, "\n")
    with_monitor(peakRAM_testfun(nsubj, ntax, d, include_ttt))
}

## start at 2.7 Mb
gc()
if (FALSE) debug(peakRAM_testfun)
## ThetaHat <- F$laplace(r)$newton(p) - up to 16 Gb, drops back to 7.3 Gb
r0 <- do_test(100)
## grows to 15.8 Mb

if (do_slow_tests) {
    res2 <- list()
    for (i in seq(nrow(testframe2))) {
        res2[[i]] <- do_test(i)
        saveRDS(res2, "reproducible/memusg2.rds")
    }
}

## still something funky with with(), I think ...
t100 <- testframe2[100,] 
attach(t100)
on.exit(detach(t100))
with_monitor(peakRAM_testfun(nsubj, ntax, d, include_ttt, use_rr = FALSE))
