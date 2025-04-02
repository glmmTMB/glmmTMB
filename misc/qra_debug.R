library(qra)
library(glmmTMB)
HCbin.cll <- readRDS("HCbin.cll.rds")
qra::extractLT(p=0.99, obj=HCbin.cll, 
               a=1:8, b=9:16, eps=0, df.t=NULL)[,-2]
