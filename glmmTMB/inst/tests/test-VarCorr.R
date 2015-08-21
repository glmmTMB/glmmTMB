require(glmmTMB)

data(Orthodont, package="nlme")
fm1 <- glmmTMB(distance ~ age + (age|Subject), data = Orthodont)
(vc <- VarCorr(fm1))  ## default print method: standard dev and corr
## both variance and std.dev.
print(vc,comp=c("Variance","Std.Dev."),digits=2)
## variance only
print(vc,comp=c("Variance"))
as.data.frame(vc)
as.data.frame(vc,order="lower.tri")
