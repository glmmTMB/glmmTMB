library(glmmTMB)
library(lme4)
data(randu)
randu$A <- factor(rep(c(1,2), 200))
randu$B <- factor(rep(c(1,2,3,4), 100))
randu$y2 <- simulate(~ x + z + (0 +x|A) + (1|B),
                     newdata = randu,
                     seed = 101,
                     newparams = list(beta = rep(1,3),
                                      theta = c(1,1),
                                      sigma = 1),
                     family = "gaussian")[[1]]

# Model
test0 <- glmmTMB(y ~ x + z + (0 +x|A) + (1|B), family="gaussian", data=randu)
test1 <- update(test0,
                start = list(theta = c(0,log(1e3))),
                map = list(theta = factor(c(1,NA))))
test2 <- update(test0,
                start = list(beta = c(1,0,0)),
                map = list(beta = factor(c(1,NA,2))))
test3 <- update(test1, y2 ~ .)
test4 <- update(test2, y2 ~ .)

## FIXME: set up metadata by crossing()/expand.grid first then iterate over rows?
results <- list()
i <- 1
for (parm in c("NULL", "beta_", "theta_")) {
    pnm <- parm
    if (parm == "NULL") parm <- NULL
    for (model in 1:4) {
        for (include_mapped in c(FALSE, TRUE)) {
            for (method in c("wald", "profile", "uniroot")) {
                cat(pnm, model, include_mapped, method, "\n")
                cc <- try(confint(get(paste0("test", model)),
                                  parm = parm,
                                  method = method,
                                  include_mapped = include_mapped), silent = TRUE)
                results[[i]] <- data.frame(parm = pnm,
                                       model = model,
                                       include_mapped = include_mapped,
                                       method = method)
                results[[i]]$status <- if (!inherits(cc, "try-error")) { "OK"
                               } else {
                                   attr(cc, "condition")$message
                               }
                i <- i + 1
                ## print(cc)
            }
        }
    }
}
results <- do.call(rbind, results)
View(results[results$status != "OK",])
## current:
## wald always OK
## beta_ OK except for uniroot + include_mapped (all models)
## {NULL, theta_} + profile works only for model 4/FALSE
##     works with 4/no-mapped
##     fails with 4/mapped (subscript out of bounds)
##     fails for other reasons for models 1, 2, 3
## {NULL, theta_} + uniroot fails ('length of dimnames [2] not equal to array extent')
##   as long as include_mapped is TRUE (all models)
