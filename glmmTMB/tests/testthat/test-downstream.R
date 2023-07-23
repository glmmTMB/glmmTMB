stopifnot(require("testthat"),
          require("glmmTMB"))

data(sleepstudy,package="lme4")
## m <- load(system.file("test_data","models.rda",package="glmmTMB", mustWork=TRUE))

## from GH 920
set.seed(101)
spcount <- rpois(120, rnorm(120, mean = rep(c(1.8, 2, 2.3), each = 40),
                            sd = 0.4))
spdat <- data.frame(species = factor(rep(1:3, each = 40)),
                  pop = factor(rep(c(1, 2, 11, 13, 21, 24), each = 20)),
                  count = spcount)

nested <- glmmTMB(count ~ species/pop, family = poisson,
                  data=spdat, control = glmmTMBControl(rank_check = "adjust"))
## equivalent glm() fit for comparison
nested0 <- glm(count ~ species + species/pop, family = poisson, data=spdat)

if (require(emmeans)) {
  test_that("emmeans", {
    skip_on_cran()
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

    ## GH780
    test3 <- glmmTMB(count ~ spp,
                     zi=~spp + mined,
                     family=truncated_nbinom2, data=Salamanders)

    e1 <- emmeans(test3,c("spp","mined"),component="zi",type="response")
    expect_is(e1, "emmGrid")
  }
  ) ## test_that

  test_that("joint_test with non-estimable terms", {
      
      expect_equal(
          c(suppressMessages(joint_tests(nested0))),
          c(suppressMessages(joint_tests(nested))),
                   tolerance = 1e-5)

      ## joint_tests(nested, component = "cmean")
  }
  )

  test_that("more on non-estimable terms", {
      excl <- with(Salamanders, which((mined == "yes") &
                                      (spp %in% c("PR","DES-L"))))
      m3 <- glmmTMB(count ~ spp * mined + (1|site),
                    zi=~spp * mined,
                    family=nbinom2, data=Salamanders[-excl, ])
      p3cond <- predict(emmeans(m3, ~spp|mined, comp = "cond"))
      p3resp <- predict(emmeans(m3, ~spp|mined, comp = "resp"))
      expect_equal(p3cond[1:3], c(-0.24506247, NA,  0.15354833))
      expect_equal(p3resp[1:3], c(0.03741431, NA, 0.40555130))
  })
} ## if require(emmeans)


if (require(effects)) {
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
