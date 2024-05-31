## exploring starting values, esp for AR1 fits
## good commit: 5939686 (right before Gaussian parameterization switch) [Oct 2023]
library(tidyverse); theme_set(theme_bw())
library(glmmTMB)
remotes::install_github("mccarthy-m-g/alda")
library(alda) # For data, not yet on CRAN: 

get_logSD <- function(obj, data) {
    ff <- reformulas::nobars(formula(obj))
    lm_fit <- lm(ff, data = data)
    return(log(abs(coef(lm_fit))))
}

sval_fun <- function(theta1, theta2, model = fm2) {
    m <- suppressWarnings(
        update(model, start = list(theta = c(theta1, theta2, 0)))
    )
    ## is.na(AIC) is a quick-and-dirty non-pos-def detector
    vv <- getME(m, "theta")
    return(data.frame(var1 = vv[1],
                      var2 = vv[2],
                      npd = is.na(AIC(m))))
}

theta_fit <- function(min = -2, max = 2, n = 21, model = fm2) {
    vals <- expand.grid(theta1 = seq(min, max, length = n),
                        theta2 = seq(min, max, length = n))
    ssv <- purrr::safely(sval_fun,
                         otherwise = data.frame(var1  = NA,
                                                var2 = NA,
                                                npd = NA))
    sv <- function(x,y) {
        ssv(x,y, model = model)$result
    }
    res <- purrr:::map2_dfr(vals$theta1, vals$theta2, sv,
                            .progress=TRUE)
    res2 <- bind_cols(vals, res)
    res2L <- pivot_longer(res2, -c(theta1, theta2))
    return(res2L)
}

theta_plots <- function(res2L, minval = -3) {
    gg1 <- ggplot(filter(res2L, name != "npd"),
           ## use 1e-10 (+log10 trans)
           ## to effectively put a floor on very small values
           aes(theta1, theta2, fill = pmax(value, minval))) +
        geom_raster() +
        facet_wrap(~name, scale = "free") +
        scale_fill_viridis_c()
    
    gg2 <- ggplot(filter(res2L, name == "npd"),
                  aes(theta1, theta2, fill = factor(value))) +
        geom_raster() +
        scale_fill_manual(values = c("black", "white"))

    list(gg1, gg2)
}

do_slow <- FALSE

##

data(sleepstudy, cbpp, Pastes,
     package = "lme4")
fsleepstudy <- transform(sleepstudy,fDays=cut(Days,c(0,3,6,10),right=FALSE),
                         row=factor(seq(nrow(sleepstudy))))

fm_ar1 <- glmmTMB(Reaction ~ 1 +
                      (1|Subject) + ar1(row+0| Subject), fsleepstudy)
VarCorr(fm_ar1)

res2L <- theta_fit(model = fm_ar1)
saveRDS(res2L, file = "fm_ar1_mstart_oldcode.rds")

theta_plots(res2L)[[1]]
theta_plots(res2L)[[2]]

v <- get_logSD(fm_ar1, fsleepstudy)
fm_ar1B <- update(fm_ar1, start = list(theta = c(v, v, 0)))

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
