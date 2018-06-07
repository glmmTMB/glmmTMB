require(glmmTMB)
require(testthat)

data(sleepstudy,package="lme4")

if (require(emmeans)) {
    context("emmeans")
    m1 <- glmmTMB(SiblingNegotiation ~ FoodTreatment*SexParent +
                (1|Nest)+offset(log(BroodSize)),
                family = nbinom1(), zi = ~1, data=Owls)
    em1 <- emmeans(m1, poly ~ FoodTreatment | SexParent)
    em2 <- emmeans(m1, poly ~ FoodTreatment | SexParent, type = "response")
    expect_is(em1,"emm_list")
    expect_true(any(grepl("given on the log (not the response) scale",
                    capture.output(print(em1)),fixed=TRUE)))
    expect_true(any(grepl("back-transformed from the log scale",
                    capture.output(print(em2)))))
    ## FIXME: need a real comparison here!
}

if (require(car)) {
    context("car::Anova")
    fm1 <- glmmTMB(Reaction~Days+(1|Subject),sleepstudy)
    ## lme4 is imported so we don't need to explicitly require() it
    fm0 <- lme4::lmer(Reaction~Days+(1|Subject),sleepstudy,REML=FALSE)
    expect_equal(Anova(fm1),Anova(fm0),tolerance=3e-6)
    expect_equal(Anova(fm1,type="III"),Anova(fm0,type="III"),tolerance=3e-6)
}

if (require(effects)) {
    context("effects")
    ## pass dd: some kind of scoping issue in testthat context
    f <- function(x,dd) {
        sapply(allEffects(x),
               function(y) {
            y$transformation$inverse(y$fit)
        })
    }
    fm2_tmb <- glmmTMB(round(Reaction)~Days+(1|Subject),family=poisson,data=sleepstudy)
    fm2_lmer <- lme4::glmer(round(Reaction)~Days+(1|Subject),family=poisson,data=sleepstudy)
    expect_equal(f(fm2_tmb),f(fm2_lmer),tolerance=2e-5)
    ## 
    set.seed(101)
    dd <<- data.frame(y=rnbinom(1000,mu=4,size=1),
                     x = rnorm(1000),
                     f=factor(rep(LETTERS[1:20],each=50)))
    fm3_tmb <- glmmTMB(y~x,family=nbinom2,data=dd)
    fm3_MASS <- MASS::glm.nb(y~x,data=dd)
    expect_equal(f(fm3_tmb,dd),f(fm3_MASS,dd),tolerance=2e-5)
}
