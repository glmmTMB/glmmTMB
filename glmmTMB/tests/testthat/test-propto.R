stopifnot(require("testthat"),
          require("glmmTMB"))

## make sure tests don't run in parallel except where we want them to
op <- getOption("glmmTMB.cores", 1)
on.exit(options(glmmTMB.cores = op))
options(glmmTMB.cores = 1)

if (require("ade4", quietly = TRUE) && require("ape", quietly = TRUE)) {
  
  set.seed(1001)
  nsp <- 50
  beta <- 1.5
  dat <- data.frame(spp = 1:nsp,
                    x = rbinom(nsp, 1, 0.4),
                    spp_int = rnorm(nsp, mean=1, sd=0.6), 
                    spp_slope = rnorm(nsp, mean=1, sd=0.6))
  dat$y <- with(dat, rpois(nsp, exp(x*beta + spp_int + x*spp_slope)))
  
  test_that("test propto fit - Identity", {
    #fit using map
    ft2 <- glmmTMB(y ~ x + diag( x | spp),  data = dat, family = poisson,
                   map = list(theta = factor(c(rep(1, 2)))))
    mat.I2 <- diag(2)
    colnames(mat.I2) <- rownames(mat.I2) <- ft2$modelInfo$reTrms$cond$cnms$spp
    ft2.I <- glmmTMB(y ~ x + propto(x | spp, mat.I2), data = dat, family = poisson)
    ## check that var matrix is as the same
    expect_equal(c(VarCorr(ft2.I)$cond[[1]]),
                 c(VarCorr(ft2)$cond[[1]]),
                 tolerance = 1e-6)
  })
  
  if (require("nlme", quietly = TRUE)) {
    # example taken from fixcorr
    data(lizards, package = "ade4")
    tree <- read.tree(text=lizards$hprA)
    liz <- lizards$traits[tree$tip.label, ]
    liz$spp <- factor(rownames(liz), levels=rownames(liz))
    liz$dummy <- factor(0)
    mat <- vcv(tree, corr=TRUE)# construct matrix
    fit <- gls(matur.L ~ age.mat,
               correlation = corSymm(mat[lower.tri(mat)], fixed=TRUE),
               data = liz)
  }
  
  test_that("test propto fit - with corr", {
    fit_phylo <- glmmTMB(matur.L ~ age.mat + propto(0 + spp | dummy, mat),
                         data = liz, dispformula = ~0)
    
    ## check equal coefficients to gls() fit
    all.equal(c(fixef(fit_phylo)$cond), c(coef(fit)),
              tolerance=5e-6)
    ## check equal coefficients to gls() fit
    all.equal(c(fixef(fit_phylo)$cond),c(coef(fit)),
              tolerance=5e-6)
    ## check that the correlation matrix is the same as constructed
    cc <- attr(VarCorr(fit_phylo)$cond[[1]], "correlation")
    dimnames(cc) = lapply(dimnames(cc), function(x) gsub("^spp","",x))
    expect_equal(cc, mat,
                 tolerance = 1e-6)
    
    })
  
    test_that("propto error about non-matrix", {
        junk <- "junk"
        expect_error( glmmTMB(y ~ x + propto(x | spp, junk), data = dat, family = poisson),
                     "expecting a matrix for propto")
    })
    ## test dimensions of matrix
    test_that("propto error with incorrect dimensions", {
      smallmat <- mat[1:10, 1:10]
      expect_error( glmmTMB(matur.L ~ age.mat + propto(0 + spp | dummy, smallmat),
                            data = liz),
                    "matrix is not the correct dimension")
    })
    ## test names of matrix
    test_that("propto error with incorrect names", {
      mattest <- mat
      colnames(mattest)[1] <- "wrong"
      rownames(mattest) <- NULL
      expect_error( glmmTMB(matur.L ~ age.mat + propto(0 + spp | dummy, mattest),
              data = liz),
              regexp = "column or row names of the matrix do not match the terms. Expecting names:.sppSa..sppSh..sppTl..sppMc..sppMy..sppPh..sppPg..sppPa..sppPb..sppPm..sppAe..sppTt..sppTs..sppZo..sppZv..sppLa..sppLs..sppLv.",
      )
    })

    test_that("propto error with unsorted names", {
        mattest <- mat
        colnames(mattest) <- rev(colnames(mattest))
        rownames(mattest) <- NULL
        expect_error(glmmTMB(matur.L ~ age.mat + propto(0 + spp | dummy, mattest),
                             data = liz),
                     "in a different order")
    })
    
    ## FIXME: test, remove if unnecessary
    options(glmmTMB.control = op) ## just in case on.exit() is inappropriate?
} ## if require("ade4") && require("ape")
