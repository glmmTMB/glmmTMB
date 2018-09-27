stopifnot(require("testthat"),
          require("glmmTMB"),
          require("MASS"))

context("weight")


set.seed(1) 
nrep = 20
nsim = 5
sdi = .1
sdii = .2
rho = -.1
slope = .8 
ni=100

dat = expand.grid(i=1:ni, rep=1:nrep , x=c(0 ,.2, .4))
RE = MASS::mvrnorm(n = ni, mu =c(0, 0), 
		Sigma = matrix(c(sdi*sdi, rho*sdi*sdii, rho*sdi*sdii ,sdii*sdii),2,2)) 
inddat = transform(dat, y=rpois(n=nrow(dat), 
                                lambda = exp(RE[i,1] + x*(slope + RE[i,2]))))

## aggdat = ddply(inddat, ~i+x+y, summarize, freq=length(rep))

aggdat = with(inddat,as.data.frame(table(i,x,y),
                                   stringsAsFactors=FALSE))
aggdat = aggdat[with(aggdat,order(i,x,y)),] ## cosmetic/match previous
aggdat = subset(aggdat,Freq>0)              ## drop zero categories
aggdat = transform(aggdat,
                   i=as.integer(i),
                   x=as.numeric(x),
                   y=as.numeric(y))
## only difference from previous is name of weights arg (Freq vs freq)

test_that("Weights can be an argument", {
    wei_glmmtmb <<- glmmTMB(y ~ x+(x|i), data=aggdat, weight=Freq,
                            family="poisson")
        expect_equal(unname(fixef(wei_glmmtmb)$cond),
                            c(-0.00907013282660578, 0.944062427131668),
                     tolerance=1e-6)
})

ind_glmmtmb <<- glmmTMB(y ~ x+(x|i), data=inddat, family="poisson")

test_that("Estimates are the same", {
	expect_equal(summary(wei_glmmtmb)$coefficients$cond, summary(ind_glmmtmb)$coefficients$cond, tolerance=1e-6)
	expect_equal(ranef(wei_glmmtmb), ranef(ind_glmmtmb), tolerance=1e-5)
	expect_equal(AIC(wei_glmmtmb), AIC(ind_glmmtmb), tolerance=1e-5)
})
