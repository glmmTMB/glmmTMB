## exploring starting values, esp for AR1 fits
## good commit: 5939686 (right before Gaussian parameterization switch) [Oct 2023]
library(tidyverse); theme_set(theme_bw())
library(furrr)
library(future)
library(progressr)
library(patchwork)
library(nlme)

## not working?
## future::plan(multisession, workers = 3)
future::plan(sequential)
library(glmmTMB)
remotes::install_github("mccarthy-m-g/alda")
library(alda) # For data, not yet on CRAN: 

get_logSD <- function(obj, data) {
    ff <- reformulas::nobars(formula(obj))
    lm_fit <- lm(ff, data = data)
    return(log(abs(coef(lm_fit))))
}

sval_fun <- function(theta1, theta2, model = fm2, fake = FALSE) {
    if (fake) return(
                  data.frame(var1 = NA, var2 = NA,
                             var3 = NA, nll = NA, npd = NA))
    m <- suppressWarnings(
        update(model, start = list(theta = c(theta1, theta2, 0)))
    )
    ## is.na(AIC) is a quick-and-dirty non-pos-def detector
    vv <- getME(m, "theta")
    return(data.frame(var1 = vv[1],
                      var2 = vv[2],
                      var3 = log(sigma(m)),
                      nll = -1*c(logLik(m)),
                      npd = is.na(AIC(m))))
}

theta_fit <- function(min = -2, max = 2, n = 21, model = fm2) {
    vals <- expand.grid(theta1 = seq(min, max, length = n),
                        theta2 = seq(min, max, length = n))
    ssv <- purrr::safely(sval_fun,
                         otherwise = sval_fun(fake=TRUE))
    sv <- function(x,y) {
        ssv(x,y, model = model)$result
    }
    res <- furrr:::future_map2_dfr(vals$theta1, vals$theta2, sv,
                            .progress=TRUE)
    res2 <- bind_cols(vals, res)
    res2L <- pivot_longer(res2, -c(theta1, theta2))
    return(res2L)
}

theta_plots <- function(res2L, minval = -3) {
    gg_logvars <- ggplot(filter(res2L, startsWith(name, "var")),
           ## use 1e-10 (+log10 trans)
           ## to effectively put a floor on very small values
           aes(theta1, theta2, fill = pmax(value, minval))) +
        geom_raster() +
        facet_wrap(~name, scale = "free") +
        scale_fill_viridis_c()

    gg_NLL <- ggplot(filter(res2L, name == "nll"),
                  aes(theta1, theta2, fill = value)) +
        geom_raster() +
        scale_fill_viridis_c()
    
    gg_npd <- ggplot(filter(res2L, name == "npd"),
                  aes(theta1, theta2, fill = factor(value))) +
        geom_raster() +
        scale_fill_manual(values = c("black", "white"))

    tibble::lst(gg_logvars, gg_NLL, gg_npd)
}

do_slow <- FALSE

##

data(sleepstudy, cbpp, Pastes,
     package = "lme4")
fsleepstudy <- transform(sleepstudy,fDays=cut(Days,c(0,3,6,10),right=FALSE),
                         row=factor(seq(nrow(sleepstudy))))

m0 <- lm(Reaction ~ 1, fsleepstudy)
sigma(m0)
## non-pos-def!
fm_ar1 <- glmmTMB(Reaction ~ 1 +
                      (1|Subject) + ar1(row+0| Subject), fsleepstudy)
VarCorr(fm_ar1)

fm_ar1B <- lme(Reaction ~ 1, random = ~1|Subject, correlation = corAR1(form = ~1|Subject),
               fsleepstudy, method = "ML")
intervals(fm_ar1B)
VarCorr(fm_ar1)

## what about a less painful example? Simulate something ...

fsleepstudy$sim <- simulate_new(~ 1 + (1|Subject) + ar1(row+0| Subject),
                                newdata=fsleepstudy,
                                newparams = list(beta=0, betad = 1, theta = c(1, 1, 1)),
                                family = gaussian,
                                seed = 101)[[1]]
fm_ar2 <- glmmTMB(sim ~ 1 +
                      (1|Subject) + ar1(row+0| Subject), fsleepstudy)



res2L <- theta_fit(model = fm_ar1)
saveRDS(res2L, file = "fm_ar1_mstart_oldcode.rds")

res2L <- readRDS("fm_ar1_mstart.rds")
filter(res2L, name == "var1" & abs(value + 0.1) < 0.1)
res2L <- readRDS("fm_ar1_mstart_oldcode.rds")
th <- theta_plots(res2L)
th[[1]] / (th[[2]] + th[[3]])

(res2L
    |> filter(name == "nll")
    |> filter(value < 900)
    |> pull(value) |> ecdf() |> plot()
)

v <- get_logSD(fm_ar1, fsleepstudy)
fm_ar1B <- update(fm_ar1, start = list(theta = c(v, v, 0)))
fm_ar1C <- update(fm_ar1, start = list(theta = c(-0.5, 4, 0), betad=4))
fm_ar1C <- update(fm_ar1, start = list(theta = c(-0.1, 4, 0), betad=4))

## https://github.com/glmmTMB/glmmTMB/issues/1036
fm2 <- glmmTMB(
    opposites_naming_score ~
        time * I(baseline_cognitive_score - 113.4571) + (time | id),
    data = opposites_naming,
    REML = TRUE
)
fm2B <- update(fm2, start = list(theta = c(2, 2, 0)))

sval_fun(1,1)

if (do_slow) {
    vals <- expand.grid(theta1 = seq(-2, 2, length = 21),
                        theta2 = seq(-2, 2, length = 21))
    res <- purrr:::map2_dfr(vals$theta1, vals$theta2, sval_fun,
                            .progress=TRUE)
    res2 <- bind_cols(vals, res)
    res2L <- pivot_longer(res2, -c(theta1, theta2))
    ggplot(filter(res2L, name != "npd"),
           ## use 1e-10 (+log10 trans)
           ## to effectively put a floor on very small values
           aes(theta1, theta2, fill = value + 1e-10)) +
        geom_raster() +
        facet_wrap(~name, scale = "free") +
        scale_fill_viridis_c(trans = "log10")
    
    ggplot(filter(res2L, name == "npd"),
           aes(theta1, theta2, fill = factor(value))) +
        geom_raster() +
        scale_fill_manual(values = c("black", "white"))
}

## a strategy for establishing starting values (most needed when
## we're using an identity link?)

clm <- get_logSD(fm2, opposites_naming)
fm2C <- update(fm2, start = list(theta = with(as.list(clm),
                                              c(`(Intercept)`, time, 0))))

stopifnot(all.equal(VarCorr(fm2B), VarCorr(fm2C), tol = 1e-4))
