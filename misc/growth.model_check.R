# -------------------- #
# growth model example
# -------------------- #

library("withr")

do_installs <- FALSE

if (do_installs) {

  ## install various combinations of old/new TMB/glmmTMB
  ## (bisection)
  inst_pkgs <- function(dir, TMB_version, glmmTMB_version) {
    if (!dir.exists(dir)) dir.create(dir, showWarnings = FALSE)
    with_libpaths(new = dir,
                  code = remotes::install_version("TMB", version = TMB_version, repos = "http://cran.us.r-project.org"))
    with_libpaths(new = dir,
                  code = remotes::install_version("glmmTMB", version = glmmTMB_version, repos = "http://cran.us.r-project.org"))
  }

  inst_pkgs("./dev_lib", "1.9.11", "1.1.9")
  inst_pkgs("./dev_lib_update", "1.9.17", "1.1.11")
  inst_pkgs("./dev_lib_intermediate", "1.9.17", "1.1.10")
  dir.create("./dev_lib_devel")
  remotes::install_github("glmmTMB/glmmTMB/glmmTMB", lib = "dev_lib_devel")
}

get_vals <- function(which = c("residuals", "model"), up2date = TRUE) {
  which <- match.arg(which)
  library(TMB)
  library(glmmTMB)
  on.exit(detach("package:glmmTMB", unload = TRUE))
  on.exit(detach("package:TMB", unload = TRUE), add = TRUE)
  mod <- get(load("./cod.growth_oldTMB.RData"))
  if (up2date) mod <- up2date(mod)
  if (which == "model") return(mod)
  residuals(mod)
}


# calculate residuals from the 'old' model with the'old' libraries
res1 = with_libpaths(new = "./dev_lib/",
                     code = get_vals())
plot(res1)

res2 = with_libpaths(new = "./dev_lib_update/",
                     code = get_vals())
max(res2)
plot(res2)

res3 = with_libpaths(new = "./dev_lib_devel/",
                     code = get_vals())
max(res3)

## MINIMAL

## hand-compute residuals
get_pars <- function(mod) with(mod$obj$env, parList(par = last.par.best))
my_resids <- function(mod, do_pkg = TRUE) {
  if (do_pkg) {
    library(TMB)
    library(glmmTMB)
    on.exit(detach("package:glmmTMB", unload = TRUE))
    on.exit(detach("package:TMB", unload = TRUE), add = TRUE)
  }
  X <- glmmTMB::getME(mod, "X")
  Z <- glmmTMB::getME(mod, "Z")
  pp <- get_pars(mod)
  pred <- with(pp, X %*% beta + Z %*% b)
  drop(model.response(model.frame(mod0)) - pred)
}

library(Matrix)
library(TMB)
devtools::load_all("~/R/pkgs/glmmTMB/glmmTMB")
mod0 <- get(load("./cod.growth_oldTMB.RData"))
max(r0 <- my_resids(mod0))
## debug(glmmTMB::up2date)
max(r1 <- my_resids(mod1 <- up2date(mod0)), do_pkg = FALSE)
max(r2 <- my_resids(mod2 <- up2date(mod0, update_gauss_disp = TRUE)))
mf <- model.frame(mod0)
weird <- which(abs(r0-r1)>0.2)  ## rows 233-236; what are these?
weird_ind <- (seq(nrow(mf)) %in% weird)+1

## residuals() gives the same answers as my_resids()
## parameters mangled somehow?

res_mat <- cbind(orig = r0,
                 no_gauss = r1,
                 no_gauss_res = residuals(mod1),
                 gauss = r2,
                 gauss_res = residuals(mod2))

pfun <- function(X, col = 1) {
  pairs(X, gap = 0,
        panel = function(x, y, ...) {
          points(x,y , col = col); abline(a=0, b=1, col = 2)
        })
}

pfun(res_mat, col = c(1,4)[weird_ind])

b_mat <-  cbind(orig = get_pars(mod0)$b,
                no_gauss = get_pars(mod1)$b,
                gauss = get_pars(mod2)$b
                )

get_last_b <- function(mod) tail(get_pars(mod)$b, 1)
get_last_b(mod0)
get_last_b(mod1)
get_last_b(mod2) ## halved (makes sense that this should be a result of var-to-SD conversion?)

pfun(b_mat)

mf[weird,]

pairs(mf, gap = 0, col = weird_ind, pch = c(1,16)[weird_ind],
      cex = weird_ind)

pairs(mf[c("cohortf", "yearf")], gap = 0, col = weird_ind, pch = c(1,16)[weird_ind],
      cex = weird_ind)

## last cohort/year?

get_pars(mod1)

mod1 = with_libpaths(new = "./dev_lib/",
                     code = get_vals(which = "model"))

mod2 = with_libpaths(new = "./dev_lib_update/",
                     code = get_vals(which = "model"))
max(my_resids(mod2))

library(TMB) ## need TMB loaded in order to return model ...
packageVersion("TMB")

## something up with Gaussian-variance updater?


## dates: 1.1.9 2024-03-20, 1.1.11 2025-04-02

# calculate residuals from the 'old' model with the 'new' libraries with 'up2date'

mod2 = with_libpaths(new = "./dev_lib_update/",
                     code = get_vals(which = "model"))
glmmTMB::up2date(mod2)



## is this a last.par vs last.par.best issue ... ???
## Z matrix is indicator-matrix only, small diffs in b shouldn't matter ...
lp0 <- mod0$obj$env$last.par
lpb0 <- mod0$obj$env$last.par.best
identical(lp0, lpb0)
all.equal(lp0, lpb0) ## small relative difference but ... ?
sort(lp0-lpb0) ## beta values are identical, b's are different

library(glmmTMB)
X <- getME(mod0, "X")
Z <- getME(mod0, "Z")
pp <- mod0$obj$env$parList()
tail(sort(res1))
tail(sort(res2))
     
res3 <- model.response(model.frame(mod0)) - drop(with(pp, X %*% beta + Z %*% b))
tail(sort(res3))  ## not identical to res1, but close

update.mod = up2date(mod1)

# compare together
par(mfrow = c(1,2),mar = c(4,4,2,2),oma = c(0,0,1,0))
plot(res2,col = "orange2",pch = 16,ylab = "residuals")
abline(h = 0,lty = 2)
points(res1,col = "lightblue4",pch = 16)
legend("topleft",legend = c("up2date","original"),col = c("orange2","lightblue4"),
       pch = 16,bty = "n")
plot(res2,res1,pch = 16,col = "lightblue3")
abline(0,1,lty = 2)
mtext("Comparison residuals",side = 3,outer = T,font = 2,line = -1,cex = 1.5)

###

library(Matrix)
##devtools::load_all("~/R/pkgs/glmmTMB/glmmTMB")
library(glmmTMB)
mod0 <- get(load("./cod.growth_oldTMB.RData"))

get_pars <- function(mod) with(mod$obj$env, parList(par = last.par.best))
my_resids <- function(mod) {
  X <- glmmTMB::getME(mod, "X")
  Z <- glmmTMB::getME(mod, "Z")
  pp <- get_pars(mod)
  pred <- with(pp, X %*% beta + Z %*% b)
  drop(model.response(model.frame(mod0)) - pred)
}

cmpfun <- function(mod) {
  ee <- mod$obj$env
  p1 <- ee$last.par.best
  p2 <- ee$parameters
  L1 <- lengths(split(p1, names(p1)))
  L2 <- lengths(p2)
  pL <- ee$parList(par = p1)
  L3 <- lengths(pL)
  list(L1, L2[L2>0], L3[L3>0],
       pL$betadisp, tail(pL$b, 1))
}

## already different orders ...
cmpfun(mod0)
mod0 <- get(load("./cod.growth_oldTMB.RData"))
cmpfun(mod0)  ## b, beta switched but doesn't seem to matter?
## [[1]]
##     b  beta betad theta 
##   121     8     1     2 
## [[2]]
##  beta     b betad theta 
##     8   121     1     2 
## [[3]]
##  beta     b betad theta 
##     8   121     1     2 
## [[4]]
## NULL
## [[5]]
## [1] -0.02734755

max(r0 <- my_resids(mod0))
## debug(glmmTMB::up2date)
mod0 <- get(load("./cod.growth_oldTMB.RData"))
max(r1 <- my_resids(mod1 <- up2date(mod0)), do_pkg = FALSE)
cmpfun(mod1)
## definitely messed up ... but orders seem OK?

debug(up2date)
LL <- lengths(obj$env$parameters)
LL[LL>0]
## current order:
##    beta        b    theta betadisp 
##       8      121        2        1
lengths(split(obj$env$last.par.best, names(obj$env$last.par.best)))
## already in a different order ... ??? how did that happen?
##        b     beta betadisp    theta 
##      121        8        1        2 

mod0 <- get(load("./cod.growth_oldTMB.RData"))
max(r2 <- my_resids(mod2 <- up2date(mod0, update_gauss_disp = TRUE)))
mf <- model.frame(mod0)
weird <- which(abs(r0-r1)>0.2)  ## rows 233-236; what are these?
weird_ind <- (seq(nrow(mf)) %in% weird)+1
