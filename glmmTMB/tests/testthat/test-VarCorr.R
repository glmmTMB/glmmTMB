stopifnot(require("testthat"),
          require("glmmTMB"),
          require("lme4"))

context("VarCorr")
##       ---------------

data("Orthodont", package="nlme")
fm1 <- glmmTMB(distance ~ age + (age|Subject), data = Orthodont)
fm1C <-   lmer(distance ~ age + (age|Subject), data = Orthodont,
               REML=FALSE) # to compare
gm1 <- glmmTMB(incidence/size ~ period + (1 | herd),
               weights=size,
               data = cbpp, family = binomial)
gm1C <-  glmer(incidence/size ~ period + (1 | herd),
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

vPix1B <- unlist(lapply(VarCorr(fmPix1B),c))
vPix1 <- unlist(lapply(VarCorr(fmPix1)[["cond"]],c))


## "manual"  (1 | Dog / Side) :
fmPix3 <- glmmTMB(pixel ~ day + I(day^2) + (day | Dog) + (1 | Dog) +
                      (1 | Side:Dog), data = Pixel)
vPix3 <- unlist(lapply(VarCorr(fmPix3)[["cond"]],c))

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
mydata <<- cbind(dd, test2)  ## GLOBAL

## The zeros in the 10 groups:
xx <- xtabs(~ a + (sim_1 == 0), mydata)

## FIXME: actually need to fit this!

## non-trivial dispersion model
data(sleepstudy, package="lme4")
fm3 <- glmmTMB(Reaction ~ Days +     (1|Subject),
               dispformula=~ Days, sleepstudy)
cc0 <- capture.output(print(fm3))
cc1 <- capture.output(print(summary(fm3)))
expect_true(any(grepl("Dispersion model:",cc0)))
expect_true(any(grepl("Dispersion model:",cc1)))


## ??? wrong context?
# not simulated this way, but returns right structure
test_that("weird variance structure", {
    mydata <- cbind(dd, test2)
    gm <- suppressWarnings(glmmTMB(sim_1 ~ 1+(b|a), zi = ~1+(b|a),
                                   data=mydata, family=poisson()))
    ## FIXME: when printed gives
    ##   Error: length(cnms) == (nc <- length(cor)) is not TRUE
})

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

getVCText <- function(obj,...) {
    c1 <- capture.output(print(obj,...))
    c1f <- read.fwf(textConnection(c1[-(1:3)]),header=FALSE,
                    fill=TRUE,
                    widths=c(10,12,9,6))
    return(c1f[,-(1:2)]) ## just value columns
}

##expect_equal(c1,c("", "Conditional model:",
##                  " Groups   Name        Std.Dev. Corr  ", 
##                  " Subject  (Intercept) 2.19409        ",
##                  "          age         0.21492  -0.581",
##                  " Residual             1.31004        "))
expect_equal(getVCText(vc),
             structure(list(V3 = c(2.1941, 0.21492, 1.31004),
                            V4 = c(NA, -0.581, NA)),
                       .Names = c("V3", "V4"),
                       class = "data.frame", row.names = c(NA, -3L)),
             tolerance=1e-5)

## both variance and std.dev.
c2 <- getVCText(vc,comp=c("Variance","Std.Dev."),digits=2)
## c2 <- capture.output(print(vc,comp=c("Variance","Std.Dev."),digits=2))
## expect_equal(c2,
##              c("", "Conditional model:",
##                " Groups   Name        Variance Std.Dev. Corr ", 
##                " Subject  (Intercept) 4.814    2.19          ",
##                "          age         0.046    0.21     -0.58", 
##                " Residual             1.716    1.31          "))
expect_equal(c2,
             structure(list(V3 = c(4.814, 0.046, 1.716), V4 = c(2.19, 0.21, 
1.31)), .Names = c("V3", "V4"), class = "data.frame", row.names = c(NA, 
-3L)))
## variance only
c3 <- getVCText(vc,comp=c("Variance"))
## c3 <- capture.output(print(vc,))
## expect_equal(c3,
##              c("", "Conditional model:",
##               " Groups   Name        Variance Corr  ", 
##               " Subject  (Intercept) 4.814050       ",
##               "          age         0.046192 -0.581",
##               " Residual             1.716203       "))
expect_equal(c3,structure(list(V3 = c(4.814071, 0.046192, 1.716208), V4 = c(NA, 
-0.581, NA)), .Names = c("V3", "V4"), class = "data.frame", row.names = c(NA, 
-3L)),
tolerance=1e-5)

if (FALSE) {  ## not yet ...
    as.data.frame(vc)
    as.data.frame(vc,order="lower.tri")
}

