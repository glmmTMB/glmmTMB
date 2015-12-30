stopifnot(require("testthat"),
          require("glmmTMB"),
          require("lme4"))

context("VarCorr")
##       ---------------

data("Orthodont", package="nlme")
fm1 <- glmmTMB(distance ~ age + (age|Subject), data = Orthodont)
fm1C <-    lmer(distance ~ age + (age|Subject), data = Orthodont,
               REML=FALSE) # to compare
gm1 <- glmmTMB(incidence/size ~ period + (1 | herd),
               weights=size,
               data = cbpp, family = binomial)
gm1C <- glmer(incidence/size ~ period + (1 | herd),
              weights=size,
              data = cbpp, family = binomial)

expect_equal(VarCorr(fm1)[["cond"]],unclass(VarCorr(fm1C)),
             tol=1e-3)
expect_equal(VarCorr(gm1)[["cond"]],unclass(VarCorr(gm1C)),
             tol=5e-3)
## have to take only last 4 lines
expect_equal(tail(capture.output(print(VarCorr(fm1),digits=3)),4),
             capture.output(print(VarCorr(fm1C),digits=3)))

data("Pixel", package="nlme")
## nPix <- nrow(Pixel)
fmPix1 <- glmmTMB(pixel ~ day + I(day^2) + (day | Dog) + (1 | Side/Dog),
                  data = Pixel)

fmPix1B <-   lmer(pixel ~ day + I(day^2) + (day | Dog) + (1 | Side/Dog),
                  data = Pixel)
## expect_equal(VarCorr(fmPix1)[["cond"]],
##           unclass(VarCorr(fmPix1B)))

## "manual"  (1 | Dog / Side) :
fmPix3 <- glmmTMB(pixel ~ day + I(day^2) + (day | Dog) + (1 | Dog) + (1 | Side:Dog), data = Pixel)

fmP1.r <- fmPix1$obj$env$report()
## str(fmP1.r)
## List of 4
##  $ corrzi: list()
##  $ sdzi  : list()
##  $ corr  :List of 3
##   ..$ : num [1, 1] 1
##   ..$ : num [1, 1] 1
##   ..$ : num [1:2, 1:2] 1 -0.598 -0.598 1
##  $ sd    :List of 3
##   ..$ : num 16.8
##   ..$ : num 9.44
##   ..$ : num [1:2] 24.83 1.73
## fmP1.r $ corr
vv <- VarCorr(fmPix1)

set.seed(12345)
dd <- data.frame(a=gl(10,100), b = rnorm(1000))
test2 <- suppressMessages(simulate(~1+(b|a), newdata=dd, family=poisson,
                  newparams= list(beta = c("(Intercept)" = 1),
                                  theta = c(1,1,1))))
## Zero-inflation : set all i.0 indices to 0:
i.0 <- sample(c(FALSE,TRUE), 1000, prob=c(.3,.7), replace=TRUE)
test2[i.0, 1] <- 0
mydata <- cbind(dd, test2)
expect_equal(head(mydata),
             structure(list(a = structure(c(1L, 1L, 1L, 1L, 1L, 1L),
               .Label = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"),
               class = "factor"), 
    b = c(0.585528817843856, 0.709466017509524, -0.109303314681054, 
    -0.453497173462763, 0.605887455840393, -1.81795596770373), 
    sim_1 = c(40, 38, 8, 0, 0, 0)),
          .Names = c("a", "b", "sim_1"),
          row.names = c("1", "2", "3", "4", "5", "6"), class = "data.frame"),
        tol=1e-5)

## The zeros in the 10 groups:
xx <- xtabs(~ a + (sim_1 == 0), mydata)
expect_equal(head(xx,3),
             structure(c(23L, 16L, 21L, 77L, 84L, 79L),
                       .Dim = c(3L, 2L),
                       .Dimnames = structure(list(
                       a = c("1", "2", "3"),
                       `sim_1 == 0` = c("FALSE", "TRUE")),
                       .Names = c("a", "sim_1 == 0")), class = "table"))

## non-trivial dispersion model
data(sleepstudy, package="lme4")
fm3 <- glmmTMB(Reaction ~ Days +     (1|Subject),
               dispformula=~ Days, sleepstudy)
cc0 <- capture.output(print(fm3))
cc1 <- capture.output(print(summary(fm3)))
expect_true(any(grepl("Dispersion model:",cc0)))
expect_true(any(grepl("Dispersion model:",cc1)))


# not simulated this way, but returns right structure
gm <- glmmTMB(sim_1 ~ 1+(b|a), zi = ~1+(b|a), data=mydata, family=poisson())
## eight updateCholesky() warnings .. which will suppress *unless* they are in the last iter.
if (FALSE) {
    str(gm.r <- gm$obj$env$report())
## List of 4
##  $ corrzi:List of 1
##   ..$ : num [1:2, 1:2] 1 0.929 0.929 1
##  $ sdzi  :List of 1
##   ..$ : num [1:2] 3.03e-05 1.87e-04
##  $ corr  :List of 1
##   ..$ : num [1:2, 1:2] 1 0.921 0.921 1
##  $ sd    :List of 1
##   ..$ : num [1:2] 0.779 1.575
}

vc <- VarCorr(fm1)  ## default print method: standard dev and corr

c1 <- capture.output(print(vc))
expect_equal(c1,c("", "Conditional model:",
                  " Groups   Name        Std.Dev. Corr  ", 
                  " Subject  (Intercept) 2.19409        ",
                  "          age         0.21492  -0.582", 
                  " Residual             1.31004        "))

## both variance and std.dev.
c2 <- capture.output(print(vc,comp=c("Variance","Std.Dev."),digits=2))
expect_equal(c2,
             c("", "Conditional model:",
               " Groups   Name        Variance Std.Dev. Corr ", 
               " Subject  (Intercept) 4.814    2.19          ",
               "          age         0.046    0.21     -0.58", 
               " Residual             1.716    1.31          "))

## variance only
c3 <- capture.output(print(vc,comp=c("Variance")))
expect_equal(c3,
             c("", "Conditional model:",
               " Groups   Name        Variance Corr  ", 
               " Subject  (Intercept) 4.814050       ",
               "          age         0.046192 -0.582", 
               " Residual             1.716203       "))

if (FALSE) {  ## not yet ...
    as.data.frame(vc)
    as.data.frame(vc,order="lower.tri")
}

