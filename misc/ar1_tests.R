library(nlme)

# generate data -----------------------------------------------------------
set.seed(1)
# https://www.r-bloggers.com/2020/02/generating-correlation-matrix-for-ar1-model/
ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                    (1:n - 1))
  rho^exponent
}

ndays <- 20
ngrps <- 10
true_cor <- 0.7
y <- MASS::mvrnorm(n = 1, 
                   mu = rep(0, ndays * ngrps), 
                   kronecker(diag(ngrps), ar1_cor(ndays, true_cor))) +
  rnorm(ndays * ngrps)  ## residual std error 1

data <- expand.grid(day = factor(1:ndays), group = LETTERS[1:ngrps]) |>
  transform(y = y,
            nday = as.numeric(day))

if (!dir.exists("older_lib")) {
  dir.create("older_lib")
  remotes::install_version("glmmTMB", version = "1.1.10", lib = "older_lib")
}

if (!dir.exists("old_lib")) {
  dir.create("old_lib")
  install.packages("glmmTMB", lib = "old_lib")
}

if (!dir.exists("new_lib")) {
  dir.create("new_lib")
  remotes::install_github("glmmTMB/glmmTMB/glmmTMB", lib = "new_lib")
}

fitfun <- function(lib) {
  library(glmmTMB, lib.loc = lib)
  on.exit(detach("package:glmmTMB", unload = TRUE))
  print(packageVersion("glmmTMB"))
  fit_glmm <- glmmTMB(y ~  1 + ar1(0 + day|group), data, REML = TRUE)
  ## have to print here, accessor/print methods go away when we unload ....
  vv <- VarCorr(fit_glmm)
  res <- c(ar1_sd = sqrt(vv$cond[[1]][1,1]),
           ar1_cor = cov2cor(vv$cond[[1]])[2,1],
           res_sd = sigma(fit_glmm))
  print(vv)
  print(logLik(fit_glmm))
  res <- c(res, nll = -1*c(logLik(fit_glmm)))
  invisible(res)
}


## ??? SOMETIMES ??? (something fragile/mutable?)
##   Groups   Name Std.Dev. Corr       
##   group    day1 1.21247  0.478 (ar1)
##   Residual      0.90617             

r_older <- fitfun("older_lib")
r_old <- fitfun("old_lib")
r_new <- fitfun("new_lib")

## all identical:

## Conditional model:
##  Groups   Name Std.Dev. Corr       
##  group    day1 0.74708  0.760 (ar1)
##  Residual      1.11011             


## how do I fit gls with a nugget?
fit_gls <- gls(y ~  1,
               correlation = corAR1(form = ~nday|group),
               data = data, method = "REML")
ar1_rho <-  coef(fit_gls$modelStruct$corStruct, uncon = FALSE)
fit_gls2 <- gls(y ~  1,
               correlation = corExp(form = ~nday|group),
               data = data, method = "REML")
exp_range <-  coef(fit_gls2$modelStruct$corStruct, uncon = FALSE)
stopifnot(all.equal(ar1_rho, exp(-1/exp_range)))
stopifnot(all.equal(logLik(fit_gls), logLik(fit_gls2)))
fit_gls3 <- gls(y ~  1,
               correlation = corExp(form = ~nday|group, nugget = TRUE),
               data = data, method = "REML")
exp_range3 <-  coef(fit_gls3$modelStruct$corStruct, uncon = FALSE)[["range"]]
exp(-1/exp_range3)  ## 0.7523

## translate range to rho?

summary(fit_gls)

## fit from asreml ---------------------------------------------------------
##  fit_asreml <- asreml(y ~ 1, random=~ ar1(day):group, data = data, maxiter = 30)
##  summary(fit_asreml)$varcomp
  ##>                   component std.error  z.ratio bound %ch
  ##> day:group         0.8304949 0.3117373 2.664085     P 0.6
  ##> day:group!day!cor 0.6730390 0.1691804 3.978232     U 0.4
  ##> units!R           1.0505278 0.2814746 3.732230     P 0.0`
}
