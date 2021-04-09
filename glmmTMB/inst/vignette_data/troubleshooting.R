library(glmmTMB)
dd <- read.csv("sherman.csv")
mod1 <- glmmTMB(formula = taeCPUE ~ Avg.Relief + (1 | site.name/reef.name), 
    data = dd, family = nbinom1, ziformula = ~1, dispformula = ~1)
mod2 <- update(mod1, ziformula=~0)
mod3 <- update(mod2, family=poisson)
mod4 <- update(mod3, family=tweedie)
save("mod1","mod2","mod3","mod4", "dd",file="troubleshooting.rda", version=2)
