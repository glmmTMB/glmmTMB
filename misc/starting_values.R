## exploring starting values, esp for AR1 fits

## https://github.com/glmmTMB/glmmTMB/issues/1036
##

data(sleepstudy, cbpp, Pastes,
     package = "lme4")
fsleepstudy <- transform(sleepstudy,fDays=cut(Days,c(0,3,6,10),right=FALSE),
                         row=factor(seq(nrow(sleepstudy))))

library(glmmTMB)
fm_ar1 <- glmmTMB(Reaction ~ 1 +
                      (1|Subject) + ar1(row+0| Subject), fsleepstudy)
VarCorr(fm_ar1)

remotes::install_github("mccarthy-m-g/alda")
library(alda) # For data, not yet on CRAN: 
fm2 <- glmmTMB(
    opposites_naming_score ~
        time * I(baseline_cognitive_score - 113.4571) + (time | id),
    data = opposites_naming,
    REML = TRUE
)
sval_fun <- function(theta1, theta2) {
    m <- suppressWarnings(
        update(fm2, start = list(theta = c(theta1, theta2, 0)))
    )
    ## is.na(AIC) is a quick-and-dirty non-pos-def detector
    vv <- unname(c(diag(VarCorr(m)$cond[[1]])))
    return(data.frame(int_var = vv[1],
                      time_var = vv[2],
                      npd = is.na(AIC(m))))
}
sval_fun(1,1)
vals <- expand.grid(theta1 = seq(-2, 2, length = 21),
                    theta2 = seq(-2, 2, length = 21))
library(tidyverse); theme_set(theme_bw())
library(patchwork)
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
       ## use 1e-10 (+log10 trans)
       ## to effectively put a floor on very small values
       aes(theta1, theta2, fill = factor(value))) +
    geom_raster() +
    scale_fill_manual(values = c("black", "white"))

