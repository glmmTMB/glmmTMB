library(glmmTMB)
devtools::load_all("~/R/pkgs/glmmTMB/glmmTMB")
library(emmeans)
library(car)
library(DHARMa)
library(performance)

options(contrasts = c("contr.sum","contr.poly"))

sp1dat <- data.frame(species = 1, pop = rep(seq(1, 2), each=20),   count = c(rpois(20, rnorm(1, 1.8, 0.4)), rpois(20, rnorm(1, 1.8, 0.4))))
sp2dat <- data.frame(species = 2, pop = rep(seq(11, 13), each=20), count = c(rpois(20, rnorm(1, 2, 0.4)),   rpois(20, rnorm(1, 2, 0.4)),
                                                                             rpois(20, rnorm(1, 2, 0.4))))
sp3dat <- data.frame(species = 3, pop = rep(seq(21, 24), each=20), count = c(rpois(20, rnorm(1, 2.3, 0.4)), rpois(20, rnorm(1, 2.3, 0.4)), rpois(20, rnorm(1, 2.3, 0.4)), rpois(20, rnorm(1, 2.3, 0.4))))
dat <- rbind(sp1dat , sp2dat , sp3dat)

dat$species <- as.factor(dat$species)
dat$pop <- as.factor(dat$pop)

## simplified version 
## set.seed(101)
## count <- rpois(120, rnorm(120, mean = rep(c(1.8, 2, 2.3), each = 40), sd = 0.4))
## dat <- data.frame(species = factor(rep(1:3, each = 40)),
##                   pop = factor(rep(c(1, 2, 11, 13, 21, 24), each = 20)),
##                   count)

nested <- glmmTMB(count ~ species/pop, family = poisson, data=dat, control = glmmTMBControl(rank_check = "adjust"))
## equivalent glm() fit for comparison
nested0 <- glm(count ~ species + species/pop, family = poisson, data=dat)
summary(nested)

try(Anova(nested, type=3)) ## OK
try(Anova(nested0, type = 3))
Anova(nested, type = 3, singular.ok = TRUE)
Anova(nested0, type = 3, singular.ok = TRUE)

Anova(nested, type=2)
Anova(nested0, type = 2)

drop1(nested, . ~ ., test="Chisq")
drop1(nested0, . ~ ., test="Chisq")

update(formula(nested), . ~ . - species)

X <- model.matrix(nested)
## drop intercept and species cols
Xr <- X[, -1*c(1, grep("^species[0-9]+$", colnames(X)))]
colnames(Xr) <- gsub(":", "", colnames(Xr))
nested2 <- update(nested, formula = reformulate(colnames(Xr), response="count"),
                  data = data.frame(count = dat$count, Xr))
pchisq(logLik(nested) - logLik(nested2), df = 2,
       lower.tail = FALSE)

joint_tests(nested)
## o@V component should have NA values removed somewhere along the line??

joint_tests(nested0)

testDispersion(nested)
simulateResiduals(fittedModel = nested, plot = T)
