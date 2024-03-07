require(glmmTMB)
require(testthat)

data(sleepstudy,package="lme4")
## m <- load(system.file("test_data","models.rda",package="glmmTMB", mustWork=TRUE))
if (require(car) && getRversion()>="3.6.0") {
    ## only testing on recent R: see comments
    ##  https://github.com/glmmTMB/glmmTMB/pull/547#issuecomment-580690208
    ##  https://github.com/glmmTMB/glmmTMB/issues/493#issuecomment-578569564
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
                     structure(list(Chisq = c(9.61229896219151, 43.6023716865556,
10.880414540324), Df = c(1, 1, 4), `Pr(>Chisq)` = c(0.00193278528775503,
4.02351135879168e-11, 0.0279413796944278)), class = c("anova",
"data.frame"), row.names = c("(Intercept)", "x", "f"), heading = c("Analysis of Deviance Table (Type III Wald chisquare tests)\n",
"Response: z")))

        expect_equal(Anova(m673, type=3, component="zi"),
                     structure(list(Chisq = c(1.66987565287872, 12.4629514668969,
0.888720602608306), Df = c(1, 1, 1), `Pr(>Chisq)` = c(0.196275190560059,
0.000415103515329875, 0.345824248063703)), class = c("anova",
"data.frame"), row.names = c("(Intercept)", "x", "y"), heading = c("Analysis of Deviance Table (Type III Wald chisquare tests)\n",
"Response: z")))

    }
)

}
