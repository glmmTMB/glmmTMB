library("glmmTMB")

data(sleepstudy, cbpp, Pastes,
     package = "lme4")

fm_noRE   <- glmmTMB(Reaction ~ Days , sleepstudy)
fm1   <- glmmTMB(Reaction ~ Days + (1| Subject), sleepstudy)
fm2   <- glmmTMB(Reaction ~ Days + (Days| Subject), sleepstudy)
fm2diag   <- glmmTMB(Reaction ~ Days + diag(Days| Subject), sleepstudy)
fm0   <- update(fm2, . ~ . -Days)
## binomial, numeric response
fm2Bn  <- update(fm2, as.numeric(Reaction>median(Reaction)) ~ .,
                 family=binomial)
## binomial, factor response
fm2Bf  <- update(fm2, factor(Reaction>median(Reaction)) ~ ., family=binomial)
fm2P  <- update(fm2, round(Reaction) ~ ., family=poisson)
fm2G  <- update(fm2, family=Gamma(link="log"))
fm2NB <- update(fm2P, family=nbinom2)
## for testing sigma() against base R
fm3   <- update(fm2, . ~ Days)
fm3G <-  update(fm3, family=Gamma(link="log"))
fm3NB <- update(fm3, round(Reaction) ~ ., family=nbinom2)
## z-i model
fm3ZIP <- update(fm2, round(Reaction) ~ ., family=poisson,
                 ziformula=~(1|Subject))
## separate-terms model
fm2diag2   <- update(fm2, . ~ Days + (1| Subject)+ (0+Days|Subject))

## model with two different grouping variables
fmP <- glmmTMB(strength ~ cask + (1|batch) + (1|sample), data=Pastes)

fm4   <- glmmTMB(Murder~Illiteracy+Population+Area+`HS Grad`,
                 data=as.data.frame(state.x77), REML = TRUE)

yb <- cbind(1:10,10)
ddb <- data.frame(y=I(yb))
ddb <- within(ddb, {
     w <- rowSums(yb)
     prop <- y[,1]/w
})
f1b <- glmmTMB(y ~ 1, family=binomial(), data=ddb)
f2b <- glm    (y ~ 1, family=binomial(), data=ddb)
f3b <- glmmTMB(prop ~ 1, weights=w, family=binomial(),
                   data=ddb)
f4b <- glmmTMB(y[,1]/w ~ 1, weights=w, family=binomial(),
                   data=ddb)

save.image(file="models.rda", version=2)
