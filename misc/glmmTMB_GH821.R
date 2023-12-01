library(glmmTMB)
dd <- read.csv("glmmTMB_GH821.csv")
ZIP3 <- glmmTMB(red~z.WTRSHDArea+Limestone+(1|Catch)+offset(log(Area)),
              family="poisson",
              ziformula = ~1,data=dd,
              na.action = 'na.fail')
formula(ZIP3)
performance::r2(ZIP3)
sessionInfo()
