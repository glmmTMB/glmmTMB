## FIXME: redundant with logic in tests/testthat/setup_makeex.R
if (!identical(Sys.getenv("NOT_CRAN"), "true")) {
  L <- load(system.file("test_data", "models.rda", package = "glmmTMB"))
} else {
  library(glmmTMB)
  save_image <- TRUE

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

gm0 <- glmmTMB(cbind(incidence, size-incidence) ~ 1 +      (1|herd),
               data = cbpp, family=binomial())
gm1 <- glmmTMB(cbind(incidence, size-incidence) ~ period + (1|herd),
               data = cbpp, family=binomial())

 ## covariance structures

  fsleepstudy <- transform(sleepstudy,fDays=cut(Days,c(0,3,6,10),right=FALSE),
                        row=factor(seq(nrow(sleepstudy))))

 ## two equivalent diagonal constructions
fm_diag1 <- glmmTMB(Reaction ~ Days + diag(Days| Subject), fsleepstudy)
fm_diag2 <- glmmTMB(Reaction ~ Days + ( 1  | Subject) + (0+Days | Subject),
               fsleepstudy)
fm_diag2_lmer <- lme4::lmer(Reaction ~ Days + ( 1  | Subject) + (0+Days | Subject),
               fsleepstudy, REML=FALSE)

fm_us1 <- glmmTMB(Reaction ~ Days + (Days| Subject), fsleepstudy)
fm_cs1 <- glmmTMB(Reaction ~ Days + cs(Days| Subject), fsleepstudy)
fm_us1_lmer <- lme4::lmer(Reaction ~ Days + ( Days  | Subject),
               fsleepstudy, REML=FALSE)

fm_cs2 <- glmmTMB(Reaction ~ Days + cs(fDays| Subject), fsleepstudy)

## these would be equivalent to a compound symmetry model with *homog* variance
fm_nest <- glmmTMB(Reaction ~ Days + (1| Subject/fDays), fsleepstudy)
fm_nest_lmer <- lme4::lmer(Reaction ~ Days + (1|Subject/fDays), fsleepstudy,
             REML=FALSE)

## model with ~ Days + ... gives non-pos-def Hessian
fm_ar1 <- glmmTMB(Reaction ~ 1 +
                      (1|Subject) + ar1(row+0| Subject), fsleepstudy)

if (save_image) save.image(file="models.rda", version=2)

 } ## if not on CRAN
