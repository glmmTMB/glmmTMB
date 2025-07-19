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
y <- MASS::mvrnorm(n = 1, 
                   mu = rep(0, ndays * ngrps), 
                   kronecker(diag(ngrps), ar1_cor(ndays, 0.7))) +
  rnorm(ndays * ngrps)  ## residual std error 1

data <- expand.grid(day = factor(1:ndays), group = LETTERS[1:ngrps]) |> 
  transform(y = y,
            nday = as.numeric(day),
            obs = factor(seq(nrow(data))))


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
  ## have to print here, print method goes away when we unload ....
  print(VarCorr(fit_glmm))
}

fitfun("older_lib")
fitfun("old_lib")
fitfun("new_lib")

## all identical:

## Conditional model:
##  Groups   Name Std.Dev. Corr       
##  group    day1 0.74708  0.760 (ar1)
##  Residual      1.11011             


## how do I fit gls with a nugget?
fit_gls <- gls(y ~  1,
               correlation = corAR1(form = ~nday|group),
               data = data, method = "REML")
fit_gls2 <- gls(y ~  1,
               correlation = corExp(form = ~nday|group),
               data = data, method = "REML")
exp(-1/0.656019) ## matches rho from fit_gls
all.equal(logLik(fit_gls), logLik(fit_gls2))
fit_gls3 <- gls(y ~  1,
               correlation = corExp(form = ~nday|group, nugget = TRUE),
               data = data, method = "REML")
exp(-1/3.513614)  ## 0.7523

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
