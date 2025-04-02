library(glmmTMB)
dd <- read.csv("sherman.csv")
mod1 <- glmmTMB(formula = taeCPUE ~ Avg.Relief + (1 | site.name/reef.name), 
    data = dd, family = nbinom1, ziformula = ~1, dispformula = ~1)
mod2 <- update(mod1, ziformula=~0)
mod3 <- update(mod2, family=poisson)
mod4 <- update(mod3, family=tweedie)
mod5 <- update(mod2, control=glmmTMBControl(optimizer=optim,
                                            optArgs=list(method="BFGS")))
mod_list <- list(base =  mod1, nozi = mod2, pois = mod3, tweedie = mod4, nozi_optim = mod5)
mod_sum <- lapply(mod_list, summary)
mod_diag <- lapply(mod_list, function(x) capture.output(diagnose(x)))
mod_pars <- lapply(mod_list, function(x) unlist(fixef(x)))
save("mod_sum", "mod_diag", "mod_pars", file="troubleshooting.rda", version=2)
