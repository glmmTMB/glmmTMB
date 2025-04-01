library(glmmTMB)
## dd <- readxl::read_excel("TMB_Data.xlsx",na="NA")
## dd <- dd[,c("taeCPUE","Avg.Relief","site.name","reef.name")]

mod1 <- glmmTMB(formula = taeCPUE ~ Avg.Relief + (1 | site.name/reef.name), 
                data = dd, family = nbinom1, ziformula = ~1, dispformula = ~1)
mod1_sum <- summary(mod1)

mod1_diag <- capture.output(diagnose(mod1))

## warning: non-integer counts ...
mod2 <- update(mod1, ziformula=~0)
mod2_sum <- summary(mod2)

mod2_diag <- capture.output(diagnose(mod1))

mod2_optim <- update(mod2,
                control=glmmTMBControl(optimizer=optim,
                                       optArgs=list(method="BFGS")))

## non-integer counts + 'false convergence'
mod3 <- update(mod2, family=poisson)
## non-integer counts
save("mod1_sum","mod2_sum","mod1_diag", "mod2_optim","mod3", version=2, file="sherman.rda")
