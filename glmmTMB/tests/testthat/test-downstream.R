require(glmmTMB)
require(testthat)

data(sleepstudy,package="lme4")
m <- load(system.file("test_data","models.rda",package="glmmTMB",
                      mustWork=TRUE))

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
    expect_equal(summary(em1[[2]])$estimate[1], -0.8586306, tolerance=1e-4)
    expect_equal(summary(em2[[2]])$ratio[1], 0.42374, tolerance=1e-4)
    
    m2 <- glmmTMB(count ~ spp + mined + (1|site),
                  zi=~spp + mined,
                  family=nbinom2, data=Salamanders)
    rgc <- ref_grid(m2, component = "cond")
    expect_is(rgc, "emmGrid")
    expect_equal(predict(rgc)[2], -1.574079, tolerance=1e-4)
    expect_equal(predict(rgc, type="response")[2], 0.207198, tolerance=1e-4)
    
    rgz <- ref_grid(m2, component = "zi")
    expect_is(rgz, "emmGrid")
    expect_equal(predict(rgz)[2], 2.071444, tolerance=1e-4)
    expect_equal(predict(rgz, type="response")[2], 0.88809654, tolerance=1e-4)

    ## test zeroing out non-focal variance components
    V <- vcov(m2)[["cond"]]
    v <- V["minedno","minedno"]
    V[] <- 0
    V["minedno","minedno"] <- v
    expect_equal(as.data.frame(emmeans(m2, ~mined, component="cond"))[["SE"]],
                 c(0.38902257366905, 0.177884950308125))
    expect_equal(as.data.frame(emmeans(m2, ~mined, component="cond", vcov.=V))[["SE"]],
                 c(0, 0.366598230362198))
}

if (require(car) && getRversion()>="3.6.0") {
    ## only testing on recent R: see comments
    ##  https://github.com/glmmTMB/glmmTMB/pull/547#issuecomment-580690208
    ##  https://github.com/glmmTMB/glmmTMB/issues/493#issuecomment-578569564
    context("car::Anova")
    fm1 <- glmmTMB(Reaction~Days+(1|Subject),sleepstudy)
    ## lme4 is imported so we don't need to explicitly require() it
    fm0 <- lme4::lmer(Reaction~Days+(1|Subject),sleepstudy,REML=FALSE)
    expect_equal(Anova(fm1),Anova(fm0),tolerance=3e-6)
    expect_equal(Anova(fm1,type="III"),Anova(fm0,type="III"),tolerance=3e-6)
    ## test Anova on various components
    fmd <- glmmTMB(Reaction~Days+(1|Subject),
                   disp=~I(Days>5), sleepstudy, REML=FALSE)
    ad <- Anova(fmd,component="disp")
    expect_equal(ad[1,1],18.767,tolerance=1e-5)
    expect_equal(rownames(ad), "I(Days > 5)")
    ac <- Anova(fmd,component="cond")
    expect_equal(ac[1,1], 160.1628, tolerance=1e-5)
    expect_equal(rownames(ac), "Days")
    expect_error(Anova(fmd,component="zi"), "trivial fixed effect")

    ## test that zi and cond tests are different
    a1 <- Anova(fm3ZIP)
    a2 <- Anova(fm3ZIP, component="zi")
    a3 <- Anova(fm3ZIP, type="III")
    a4 <- Anova(fm3ZIP, type="III",component="zi")
    get_pval <- function(x) c(na.omit(x$`Pr(>Chisq)`))
    expect_equal(get_pval(a1),1.82693150434104e-13)
    expect_equal(get_pval(a2),numeric(0))
    expect_equal(get_pval(a3),c(0, 1.82693150434104e-13))
    expect_equal(get_pval(a4),0.81337346580467)

    test_that("Anova matches zi attributes correctly", {
        ## zi and cond for cases with different models (GH 673)
        ## set up case where one of the 'term' indices matches up with an element
        ## in the 'assign' attribute for the conditional model > the number of terms for the zi model
        ## here conditional model
        ## ?? not sure why the data simulation/model fitting has to be within test_that() but apparently it does ??
        set.seed(101)
        n <- 100
        ## give data unique name to avoid interfering with dd below (only relevant in some contexts?)
        dd673 <<- data.frame(x=rnorm(n),y=rnorm(n),f=factor(sample(5,size=n,replace=TRUE)))
        beta <- rep(1,6)
        X <- model.matrix(~x+f,data=dd673)
        dd673 <<- within(dd673, {
            eta_cond <- exp(X %*% beta)
            eta_zi <- plogis(x+y)
            z <- ifelse(runif(n)<eta_zi,0,rnbinom(n,mu=eta_cond,size=1))
        })
        m673 <<- glmmTMB(z~x+f, family=nbinom2, data=dd673, zi=~x+y)
        ac <-  attr(model.matrix(m673), "assign")
        az <- attr(model.matrix(m673, component="zi"), "assign")

        expect_true(max(match(az,ac))>length(az)-1)
        expect_equal(Anova(m673, type=3),
                     structure(list(Chisq = c(4.72184220080362, 26.0412437493939, 
2.44974442018967), Df = c(1, 1, 4), `Pr(>Chisq)` = c(0.0297818205083957, 
3.34201109860344e-07, 0.653656869363576)), class = c("anova", 
"data.frame"), row.names = c("(Intercept)", "x", "f"), heading = c("Analysis of Deviance Table (Type III Wald chisquare tests)\n", 
"Response: z")))
        expect_equal(Anova(m673, type=3, component="zi"),
                     structure(list(Chisq = c(0.605706252119866, 4.82515004746816, 
8.93942508493268), Df = c(1, 1, 1), `Pr(>Chisq)` = c(0.436409077669975, 
0.0280474237991421, 0.00279080545380766)), class = c("anova", 
"data.frame"), row.names = c("(Intercept)", "x", "y"), heading = c("Analysis of Deviance Table (Type III Wald chisquare tests)\n", 
"Response: z")))
    }
)

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
    if (getRversion() >= "3.6.0") {
        ## only testing on recent R: see comments
        ##  https://github.com/glmmTMB/glmmTMB/pull/547#issuecomment-580690208
        ##  https://github.com/glmmTMB/glmmTMB/issues/493#issuecomment-578569564
        expect_equal(f(fm2_tmb),f(fm2_lmer),tolerance=2e-5)
        ## 
        set.seed(101)
        dd <<- data.frame(y=rnbinom(1000,mu=4,size=1),
                          x = rnorm(1000),
                          f=factor(rep(LETTERS[1:20],each=50)))
        fm3_tmb <- glmmTMB(y~x,family=nbinom2,data=dd)
        fm3_MASS <- MASS::glm.nb(y~x,data=dd)
        ## suppressing "overriding variance function for effects: computed variances may be incorrect" warning here
        expect_equal(suppressWarnings(f(fm3_tmb,dd)),f(fm3_MASS,dd),tolerance=2e-5)
    } ## recent R
} ## effects
