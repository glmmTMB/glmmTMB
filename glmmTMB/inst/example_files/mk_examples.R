library(glmmTMB)
salamander1 <- glmmTMB(count ~ mined + (1|site),
                       zi=~mined,
                       family=poisson, data=Salamanders)
saveRDS(salamander1,"salamander1.rds",version=2)
