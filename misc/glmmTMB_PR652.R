## utility for checking parallel issues
## https://github.com/glmmTMB/glmmTMB/issues/620
## https://github.com/glmmTMB/glmmTMB/pull/652

tag <- "BMB"     ## indicator of who is running/ran these examples
ref <- "master"    ## what branch/package version?
## try with "CRAN", "master", "#652"

## adjust these to your liking
max_cores <- min(16,parallel::detectCores()-1)
cores_vec <- seq(max_cores)
reps_vec <- c(1,5,10)

fn <- sprintf("glmmTMB_PR652_%s_%s",ref,tag)
datfn <- sprintf("%s.rds",fn)
chkfn <- sprintf("%s_chk.RData",fn)

## fresh package install from specified source
if (ref=="CRAN") {
    install.packages("glmmTMB",repos="https://cloud.r-project.org")
} else {
    remotes::install_github("glmmTMB/glmmTMB/glmmTMB",ref=ref)
}

if (length(fn)==1) {  ## doing a run rather than compiling results
    ## run
    data("InstEval", package = "lme4")
    library(glmmTMB)

    f1 <- function(p, r=1, seed=101, debug=TRUE) {
        if (debug) cat(sprintf("fitting %d rep(s) with %d core(s)\n",r,p))
        set.seed(seed)
        d <- InstEval[rep(seq(nrow(InstEval)),r),]
        form <- y ~ service + lectage + studage + (1 | d) + (1 | s)
        ctrl <- glmmTMBControl(parallel=p)
        tt <- system.time(res <- try(glmmTMB(form, data = d, control = ctrl),
                                     silent=TRUE))
        ## if error, save the error message
        res_str <- if (inherits(res,"try-error")) {
                       attr(res,"condition")$message
                   } else "OK"
        print(res_str)
        print(tt)
        return(list(result=res_str,time=tt))
    }

    ## loop over number of times InstEval data is replicated
    times <- vector("list", length(reps_vec))
    names(times) <- paste0("rep_",reps_vec)
    attr(times,"sessionInfo") <- sessionInfo()
    for (i in seq_along(reps_vec)) {
        times[[i]] <- vector("list", length(cores_vec))
        names(times[[i]]) <- paste0("cores_",cores_vec)
        ## loop over number of cores to use
        for (j in seq_along(cores_vec)) {
            cat(i,j,"reps:",reps_vec[[i]]," cores:",cores_vec[[j]],"\n")
            times[[i]][[j]] <- f1(p=cores_vec[[j]],r=reps_vec[[i]])
            save.image(file=chkfn)  ## checkpoint
        }
    }
    saveRDS(times,file=datfn)
}

### compilation of results ...
if (FALSE) {
    load(chkfn)
    ## OBSOLETE
    load(outfile)
    tmpf <- function(z) {
        t(sapply(z, function(x) {
            if (is(x,"character")) return(rep(NA,2))
            x[c("user.self","elapsed")]
        }))
    }
    tmpf2 <- function(z) {
        sapply(z, function(x) {
            if (!is(x,"character")) return(NA)
            return(x)
        })
    }
    
    cbind(cores=1:16,tmpf(times1))
    cbind(cores=1:16,tmpf(times2))
}


