stopifnot(require("testthat"),
          require("glmmTMB"),
          require("lme4"))

data(sleepstudy, cbpp,
     package = "lme4")
sleepstudy <- transform(sleepstudy,fDays=cut(Days,c(0,3,6,10),right=FALSE),
                        row=factor(seq(nrow(sleepstudy))))

context("variance structures")

## two equivalent diagonal constructions
fm_diag1 <- glmmTMB(Reaction ~ Days + diag(Days| Subject), sleepstudy)
fm_diag2 <- glmmTMB(Reaction ~ Days + ( 1  | Subject) + (0+Days | Subject),
               sleepstudy)
fm_diag2_lmer <- lmer(Reaction ~ Days + ( 1  | Subject) + (0+Days | Subject),
               sleepstudy, REML=FALSE)

fm_us1 <- glmmTMB(Reaction ~ Days + (Days| Subject), sleepstudy)
fm_cs1 <- glmmTMB(Reaction ~ Days + cs(Days| Subject), sleepstudy)
fm_us1_lmer <- lmer(Reaction ~ Days + ( Days  | Subject),
               sleepstudy, REML=FALSE)

fm_cs2 <- glmmTMB(Reaction ~ Days + cs(fDays| Subject), sleepstudy)

## these would be equivalent to a compound symmetry model with *homog* variance
fm_nest <- glmmTMB(Reaction ~ Days + (1| Subject/fDays), sleepstudy)
fm_nest_lmer <- lmer(Reaction ~ Days + (1|Subject/fDays), sleepstudy,
             REML=FALSE)

## model with ~ Days + ... gives non-pos-def Hessian
fm_ar1 <- glmmTMB(Reaction ~ 1 +
                   ar1(row+0| Subject), sleepstudy)

test_that("diag", {
   ## two formulations of diag and lme4 all give same log-lik
   expect_equal(logLik(fm_diag1),logLik(fm_diag2_lmer))
   expect_equal(logLik(fm_diag1),logLik(fm_diag2))
})

test_that("cs_us", {
## for a two-level factor, compound symmetry and unstructured
##  give same result
    expect_equal(logLik(fm_us1),logLik(fm_cs1))
    expect_equal(logLik(fm_us1),logLik(fm_us1_lmer))
})

test_that("cs_homog", {
## *homogenous* compound symmetry vs. nested random effects
    expect_equal(logLik(fm_nest),logLik(fm_nest_lmer))

})

test_that("ar1", {
    vv <- VarCorr(fm_ar1)[["cond"]]
    cc <- cov2cor(vv[[1]])
    expect_equal(cc[1,],cc[,1])
    expect_equal(unname(cc[1,]),
                 cc[1,2]^(0:(nrow(cc)-1)))
})

## FIXME: simpler to check formatVC() directly?
get_vcout <- function(x,g="Subject") {
    cc <- capture.output(print(VarCorr(x)))
    cc1 <- grep(g,cc,value=TRUE)
    ss <- strsplit(cc1,"[^[:alnum:][:punct:]]+")[[1]]
    return(ss[nchar(ss)>0])
}

test_that("varcorr_print", {
    ss <- get_vcout(fm_cs1)
    expect_equal(length(ss),5)
    expect_equal(ss[4:5],c("0.081","(cs)"))
    ss2 <- get_vcout(fm_ar1)
    expect_equal(length(ss2),5)
    expect_equal(ss2[4:5],c("0.853","(ar1)"))

    ## test case with two different size V-C
    set.seed(101)
    dd <- data.frame(y=rnorm(1000),c=factor(rep(1:2,500)),
                 w=factor(rep(1:10,each=100)),
                 s=factor(rep(1:10,100)))
    ## non-pos-def case (we don't care at the moment)
    m1 <- suppressWarnings(glmmTMB(y~c+(c|w)+(1|s),data=dd,
                  family=gaussian))
    cc <- capture.output(print(VarCorr(m1),digits=2))
    expect_equal(cc,
    c("", "Conditional model:", " Groups   Name        Std.Dev. Corr", 
" w        (Intercept) 3.1e-05      ", "          c2          4.9e-06  0.98", 
" s        (Intercept) 3.4e-05      ", " Residual             9.6e-01      "
))
    })

