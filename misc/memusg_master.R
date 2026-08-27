source("reproducible/leverage_funs.R")

## reinstalling ...
if (FALSE) {

    remotes::install_github("mrc-ide/memprof")
    
    pkgdir <- "~/R/pkgs/glmmTMB/glmmTMB"
    git_status <- system(sprintf("git -C %s status", pkgdir), intern = TRUE)
    if (git_status[1] != "On branch formula_env") stop("need formula_env branch")
    install_glmmTMB(pkgdir, "glmmTMB_lev")
}

library(memprof)
library(glmmTMB, lib = "glmmTMB_lev")
library(peakRAM)
library(RTMB)
library(Matrix)
library(tidyverse)
do_simple_tests <- FALSE

if (do_simple_tests) {
    ## basic tests
    peakRAM_testfun(5, 5, 2)
    peakRAM_testfun(5, 5, 3)
    peakRAM_testfun(5, 5, 3, include_ttt = TRUE)
}

testframe <- expand.grid(nsubj = seq(5, 40, by = 5),
                         ntax = seq(25, 200, by = 25))
write.csv(testframe, "reproducible/testvals.csv")

testframe2 <- expand.grid(nsubj = seq(5, 40, by = 5),
                         ntax = seq(25, 200, by = 25),
                         d = 1:3,
                         include_ttt = c(TRUE, FALSE)) |>
    mutate(ntheta = ntheta(nsubj, ntax, d, include_ttt),
           nobs = ntax*nsubj) |>
    ## order by expected mem usg, to try to get in as many as possible
    ##  before out-of-memory crash
    arrange(ntheta*nobs)
write.csv(testframe2, "reproducible/testvals2.csv")

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

res2 <- list()
## something wonky with with()/i?
## what's here is clunky but ...
for (i in seq(nrow(testframe2))) {
    nsubj <- testframe2$nsubj[i]
    ntax <-  testframe2$ntax[i]
    d <- testframe2$d[i]
    include_ttt <- testframe2$include_ttt[i]
    ## with(testframe2[i,], {
        cat(i, nsubj, ntax, d, include_ttt, "\n")
        res2[[i]] <- with_monitor(peakRAM_testfun(nsubj, ntax, d, include_ttt))
        saveRDS(res2, "reproducible/memusg2.rds")
    ##})
}
         

